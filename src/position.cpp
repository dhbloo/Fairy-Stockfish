/*
  Stockfish, a UCI chess playing engine derived from Glaurung 2.1
  Copyright (C) 2004-2022 The Stockfish developers (see AUTHORS file)

  Stockfish is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  Stockfish is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <algorithm>
#include <cassert>
#include <cstddef> // For offsetof()
#include <cstring> // For std::memset, std::memcmp
#include <iomanip>
#include <sstream>

#include "bitboard.h"
#include "misc.h"
#include "movegen.h"
#include "position.h"
#include "thread.h"
#include "tt.h"
#include "uci.h"
#include "syzygy/tbprobe.h"

using std::string;

namespace Stockfish {

    namespace Zobrist {

        Key psq[PIECE_NB][SQUARE_NB];
        Key side, noPawns;
        Key inHand[PIECE_NB][SQUARE_NB];
        Key checks[COLOR_NB][CHECKS_NB];
    }


    /// operator<<(Position) returns an ASCII representation of the position

    std::ostream& operator<<(std::ostream& os, const Position& pos) {

        os << "\n ";
        for (File f = FILE_A; f <= pos.max_file(); ++f)
            os << "+---";
        os << "+\n";

        for (Rank r = pos.max_rank(); r >= RANK_1; --r)
        {
            for (File f = FILE_A; f <= pos.max_file(); ++f)
                if (pos.unpromoted_piece_on(make_square(f, r)))
                    os << " |+" << pos.piece_to_char()[pos.unpromoted_piece_on(make_square(f, r))];
                else
                    os << " | " << pos.piece_to_char()[pos.piece_on(make_square(f, r))];

            os << " |" << (1 + r);
            if (r == pos.max_rank() || r == RANK_1)
            {
                Color c = r == RANK_1 ? WHITE : BLACK;
                if (c == pos.side_to_move())
                    os << " *";
                else
                    os << "  ";
            }
            os << "\n ";
            for (File f = FILE_A; f <= pos.max_file(); ++f)
                os << "+---";
            os << "+\n";
        }

        for (File f = FILE_A; f <= pos.max_file(); ++f)
            os << "   " << char('a' + f);
        os << "\n";
        os << "\nFen: " << pos.fen() << "\nSfen: " << pos.fen(true) << "\nKey: " << std::hex << std::uppercase
            << std::setfill('0') << std::setw(16) << pos.key()
            << std::setfill(' ') << std::dec << "\nCheckers: ";

        for (Bitboard b = pos.checkers(); b; )
            os << UCI::square(pos, pop_lsb(b)) << " ";

        os << "\nChased: ";
        for (Bitboard b = pos.state()->chased; b; )
            os << UCI::square(pos, pop_lsb(b)) << " ";

        return os;
    }

    /// Position::init() initializes at startup the various arrays used to compute hash keys

    void Position::init() {

        PRNG rng(1070372);

        for (Color c : {WHITE, BLACK})
            for (PieceType pt = PAWN; pt <= KING; ++pt)
                for (Square s = SQ_A1; s <= SQ_MAX; ++s)
                    Zobrist::psq[make_piece(c, pt)][s] = rng.rand<Key>();

        Zobrist::side = rng.rand<Key>();
        Zobrist::noPawns = rng.rand<Key>();

        for (Color c : {WHITE, BLACK})
            for (int n = 0; n < CHECKS_NB; ++n)
                Zobrist::checks[c][n] = rng.rand<Key>();

        for (Color c : {WHITE, BLACK})
            for (PieceType pt = PAWN; pt <= KING; ++pt)
                for (int n = 0; n < SQUARE_NB; ++n)
                    Zobrist::inHand[make_piece(c, pt)][n] = rng.rand<Key>();
    }


    /// Position::set() initializes the position object with the given FEN string.
    /// This function is not very robust - make sure that input FENs are correct,
    /// this is assumed to be the responsibility of the GUI.

    Position& Position::set(const Variant* v, const string& fenStr, StateInfo* si, Thread* th, bool sfen) {
        /*
           A FEN string defines a particular position using only the ASCII character set.

           A FEN string contains six fields separated by a space. The fields are:

           1) Piece placement (from white's perspective). Each rank is described, starting
              with rank 8 and ending with rank 1. Within each rank, the contents of each
              square are described from file A through file H. Following the Standard
              Algebraic Notation (SAN), each piece is identified by a single letter taken
              from the standard English names. White pieces are designated using upper-case
              letters ("PNBRQK") whilst Black uses lowercase ("pnbrqk"). Blank squares are
              noted using digits 1 through 8 (the number of blank squares), and "/"
              separates ranks.

           2) Active color. "w" means white moves next, "b" means black.

           3) Castling availability. If neither side can castle, this is "-". Otherwise,
              this has one or more letters: "K" (White can castle kingside), "Q" (White
              can castle queenside), "k" (Black can castle kingside), and/or "q" (Black
              can castle queenside).

           4) En passant target square (in algebraic notation). If there's no en passant
              target square, this is "-". If a pawn has just made a 2-square move, this
              is the position "behind" the pawn. Following X-FEN standard, this is recorded only
              if there is a pawn in position to make an en passant capture, and if there really
              is a pawn that might have advanced two squares.

           5) Halfmove clock. This is the number of halfmoves since the last pawn advance
              or capture. This is used to determine if a draw can be claimed under the
              fifty-move rule.

           6) Fullmove number. The number of the full move. It starts at 1, and is
              incremented after Black's move.
        */

        unsigned char col, row, token;
        size_t idx;
        std::istringstream ss(fenStr);

        std::memset(this, 0, sizeof(Position));
        std::memset(si, 0, sizeof(StateInfo));
        st = si;

        var = v;

        ss >> std::noskipws;

        Square sq = SQ_A1 + max_rank() * NORTH;

        // 1. Piece placement
        while ((ss >> token) && !isspace(token))
        {
            if (isdigit(token))
            {
                if (isdigit(ss.peek()))
                {
                    sq += 10 * (token - '0') * EAST;
                    ss >> token;
                }

                sq += (token - '0') * EAST; // Advance the given number of files
            }

            else if (token == '/')
            {
                sq += 2 * SOUTH + (FILE_MAX - max_file()) * EAST;
                if (!is_ok(sq))
                    break;
            }

            else if ((idx = piece_to_char().find(token)) != string::npos || (idx = piece_to_char_synonyms().find(token)) != string::npos)
            {
                if (ss.peek() == '~')
                    ss >> token;
                put_piece(Piece(idx), sq, token == '~');
                ++sq;
            }

            // Stop before pieces in hand
            else if (token == '[')
                break;
        }
        // Pieces in hand
        if (!isspace(token))
            while ((ss >> token) && !isspace(token))
            {
                if (token == ']')
                    continue;
                else if ((idx = piece_to_char().find(token)) != string::npos)
                    add_to_hand(Piece(idx));
            }

        // 2. Active color
        ss >> token;
        sideToMove = (token != (sfen ? 'w' : 'b') ? WHITE : BLACK);  // Invert colors for SFEN
        ss >> token;

        // 3-4. Skip parsing castling and en passant flags if not present
        st->epSquare = SQ_NONE;

        // Check counter for nCheck
        ss >> std::skipws >> token >> std::noskipws;
        ss.putback(token);

        // 5-6. Halfmove clock and fullmove number
        if (sfen)
        {
            // Pieces in hand for SFEN
            int handCount = 1;
            while ((ss >> token) && !isspace(token))
            {
                if (token == '-')
                    continue;
                else if (isdigit(token))
                {
                    handCount = token - '0';
                    while (isdigit(ss.peek()) && ss >> token)
                        handCount = 10 * handCount + (token - '0');
                }
                else if ((idx = piece_to_char().find(token)) != string::npos)
                {
                    for (int i = 0; i < handCount; i++)
                        add_to_hand(Piece(idx));
                    handCount = 1;
                }
            }
            // Move count is in ply for SFEN
            ss >> std::skipws >> gamePly;
            gamePly = std::max(gamePly - 1, 0);
        }
        else
        {
            ss >> std::skipws >> st->rule50 >> gamePly;

            // Convert from fullmove starting from 1 to gamePly starting from 0,
            // handle also common incorrect FEN with fullmove = 0.
            gamePly = std::max(2 * (gamePly - 1), 0) + (sideToMove == BLACK);
        }
        thisThread = th;
        set_state(st);

        assert(pos_is_ok());

        return *this;
    }


    /// Position::set_check_info() sets king attacks to detect if a move gives check

    void Position::set_check_info(StateInfo* si) const {

        si->blockersForKing[WHITE] = slider_blockers(pieces(BLACK), square<KING>(WHITE), si->pinners[BLACK], BLACK);
        si->blockersForKing[BLACK] = slider_blockers(pieces(WHITE), square<KING>(BLACK), si->pinners[WHITE], WHITE);

        Square ksq = square<KING>(~sideToMove);

        // For unused piece types, the check squares are left uninitialized
        si->nonSlidingRiders = 0;
        for (PieceType pt : pieceTypes)
        {
            si->checkSquares[pt] = ksq != SQ_NONE ? attacks_bb(~sideToMove, pt, ksq, pieces()) : Bitboard(0);
            // Collect special piece types that require slower check and evasion detection
            if (AttackRiderTypes[pt] & NON_SLIDING_RIDERS)
                si->nonSlidingRiders |= pieces(pt);
        }
        si->shak = si->checkersBB & byTypeBB[ROOK];
        si->chased = chased();
        si->legalCapture = NO_VALUE;
    }


    /// Position::set_state() computes the hash keys of the position, and other
    /// data that once computed is updated incrementally as moves are made.
    /// The function is only used when a new position is set up, and to verify
    /// the correctness of the StateInfo data when running in debug mode.

    void Position::set_state(StateInfo* si) const {

        si->key = si->materialKey = 0;
        si->pawnKey = Zobrist::noPawns;
        si->nonPawnMaterial[WHITE] = si->nonPawnMaterial[BLACK] = VALUE_ZERO;
        si->checkersBB = attackers_to(square<KING>(sideToMove), ~sideToMove);
        si->move = MOVE_NONE;

        set_check_info(si);

        for (Bitboard b = pieces(); b; )
        {
            Square s = pop_lsb(b);
            Piece pc = piece_on(s);
            si->key ^= Zobrist::psq[pc][s];

            if (type_of(pc) != KING)
                si->nonPawnMaterial[color_of(pc)] += PieceValue[MG][pc];
        }

        if (sideToMove == BLACK)
            si->key ^= Zobrist::side;

        for (Color c : {WHITE, BLACK})
            for (PieceType pt = PAWN; pt <= KING; ++pt)
            {
                Piece pc = make_piece(c, pt);

                for (int cnt = 0; cnt < pieceCount[pc]; ++cnt)
                    si->materialKey ^= Zobrist::psq[pc][cnt];
            }
    }


    /// Position::set() is an overload to initialize the position object with
    /// the given endgame code string like "KBPKN". It is mainly a helper to
    /// get the material key out of an endgame code.

    Position& Position::set(const string& code, Color c, StateInfo* si) {

        assert(code[0] == 'K');

        string sides[] = { code.substr(code.find('K', 1)),      // Weak
                           code.substr(0, std::min(code.find('v'), code.find('K', 1))) }; // Strong

        assert(sides[0].length() > 0 && sides[0].length() < 8);
        assert(sides[1].length() > 0 && sides[1].length() < 8);

        std::transform(sides[c].begin(), sides[c].end(), sides[c].begin(), tolower);

        string n = std::to_string(FILE_NB);
        string fenStr = n + "/" + sides[0] + char(FILE_NB - sides[0].length() + '0') + "/" + n + "/" + n + "/" + n + "/"
            + n + "/" + sides[1] + char(FILE_NB - sides[1].length() + '0') + "/" + n + " w - - 0 10";

        return set(variants.find("fairy")->second, fenStr, si, nullptr);
    }


    /// Position::fen() returns a FEN representation of the position. In case of
    /// Chess960 the Shredder-FEN notation is used. This is mainly a debugging function.

    string Position::fen(bool sfen, bool showPromoted, int countStarted, std::string holdings) const {

        int emptyCnt;
        std::ostringstream ss;

        for (Rank r = max_rank(); r >= RANK_1; --r)
        {
            for (File f = FILE_A; f <= max_file(); ++f)
            {
                for (emptyCnt = 0; f <= max_file() && empty(make_square(f, r)); ++f)
                    ++emptyCnt;

                if (emptyCnt)
                    ss << emptyCnt;

                if (f <= max_file())
                {
                    if (unpromoted_piece_on(make_square(f, r)))
                        // Promoted shogi pieces, e.g., +r for dragon
                        ss << "+" << piece_to_char()[unpromoted_piece_on(make_square(f, r))];
                    else
                        ss << piece_to_char()[piece_on(make_square(f, r))];
                }
            }

            if (r > RANK_1)
                ss << '/';
        }

        // SFEN
        if (sfen)
        {
            ss << (sideToMove == WHITE ? " b " : " w ");
            for (Color c : {WHITE, BLACK})
                for (PieceType pt = KING; pt >= PAWN; --pt)
                    if (pieceCountInHand[c][pt] > 0)
                    {
                        if (pieceCountInHand[c][pt] > 1)
                            ss << pieceCountInHand[c][pt];
                        ss << piece_to_char()[make_piece(c, pt)];
                    }
            if (count_in_hand(ALL_PIECES) == 0)
                ss << '-';
            ss << " " << gamePly + 1;
            return ss.str();
        }

        ss << (sideToMove == WHITE ? " w " : " b ");
        ss << '-';

        ss << (ep_square() == SQ_NONE ? " - " : " " + UCI::square(*this, ep_square()) + " ");

        ss << st->rule50;
        ss << " " << 1 + (gamePly - (sideToMove == BLACK)) / 2;

        return ss.str();
    }


    /// Position::slider_blockers() returns a bitboard of all the pieces (both colors)
    /// that are blocking attacks on the square 's' from 'sliders'. A piece blocks a
    /// slider if removing that piece from the board would result in a position where
    /// square 's' is attacked. For example, a king-attack blocking piece can be either
    /// a pinned or a discovered check piece, according if its color is the opposite
    /// or the same of the color of the slider.

    Bitboard Position::slider_blockers(Bitboard sliders, Square s, Bitboard& pinners, Color c) const {

        Bitboard blockers = 0;
        pinners = 0;

        if (s == SQ_NONE || !sliders)
            return blockers;

        // Snipers are sliders that attack 's' when a piece and other snipers are removed
        Bitboard snipers = 0;
        Bitboard slidingSnipers = 0;

        auto add_snipers = [&](PieceType pt) {
            Bitboard b = sliders & (PseudoAttacks[~c][pt][s] ^ LeaperAttacks[~c][pt][s]) & pieces(c, pt);
            if (b)
            {
                // Consider asymmetrical moves (e.g., horse)
                if (AttackRiderTypes[pt] & ASYMMETRICAL_RIDERS)
                    snipers |= b & ~attacks_by_horse(s, pieces());
                else
                    snipers |= b & ~attacks_bb(~c, pt, s, pieces());
                if (AttackRiderTypes[pt] & ~HOPPING_RIDERS)
                    slidingSnipers |= snipers & pieces(pt);
            }
        };

        add_snipers(PieceType::ROOK);
        add_snipers(PieceType::FERS);
        add_snipers(PieceType::CANNON);
        add_snipers(PieceType::SOLDIER);
        add_snipers(PieceType::HORSE);
        add_snipers(PieceType::ELEPHANT);
        add_snipers(PieceType::KING);

        Bitboard occupancy = pieces() ^ slidingSnipers;

        while (snipers)
        {
            Square sniperSq = pop_lsb(snipers);
            bool isHopper = AttackRiderTypes[type_of(piece_on(sniperSq))] & HOPPING_RIDERS;
            Bitboard b = between_bb(s, sniperSq, type_of(piece_on(sniperSq))) & (isHopper ? (pieces() ^ sniperSq) : occupancy);

            if (b && (!more_than_one(b) || (isHopper && popcount(b) == 2)))
            {
                blockers |= b;
                if (b & pieces(color_of(piece_on(s))))
                    pinners |= sniperSq;
            }
        }
        return blockers;
    }


    /// Position::attackers_to() computes a bitboard of all pieces which attack a
    /// given square. Slider attacks use the occupied bitboard to indicate occupancy.

    Bitboard Position::attackers_to(Square s, Bitboard occupied, Color c) const {

        Bitboard b = 0;

        b |= attacks_bb<ROOK>(s, occupied) & pieces(c, ROOK);
        if (board_bb(c, FERS) & s)
            b |= attacks_bb<FERS>(s, occupied) & pieces(c, FERS);
        b |= attacks_bb<CANNON>(s, occupied) & pieces(c, CANNON);
        b |= PseudoAttacks[~c][SOLDIER][s] & pieces(c, SOLDIER);
         b |= attacks_by_horse(s, occupied) & pieces(c, HORSE);
        if (board_bb(c, ELEPHANT) & s)
            b |= attacks_bb<ELEPHANT>(s, occupied) & pieces(c, ELEPHANT);
        if (board_bb(c, KING) & s)
            b |= PseudoAttacks[~c][WAZIR][s] & pieces(c, KING);

        // Unpromoted soldiers
        if (b & pieces(SOLDIER) && relative_rank(c, s, max_rank()) < RANK_6)
            b ^= b & pieces(SOLDIER) & ~PseudoAttacks[~c][SHOGI_PAWN][s];

        return b;
    }


    Bitboard Position::attackers_to(Square s, Bitboard occupied) const {
        return attackers_to(s, occupied, WHITE) | attackers_to(s, occupied, BLACK);
    }

    /// Position::attackers_to_pseudo_royals computes a bitboard of all pieces
    /// of a particular color attacking at least one opposing pseudo-royal piece
    Bitboard Position::attackers_to_pseudo_royals(Color c) const {
        Bitboard attackers = 0;
        Bitboard pseudoRoyals = st->pseudoRoyals & pieces(~c);
        Bitboard pseudoRoyalsTheirs = st->pseudoRoyals & pieces(c);
        while (pseudoRoyals) {
            Square sr = pop_lsb(pseudoRoyals);
            attackers |= attackers_to(sr, c);
        }
        return attackers;
    }


    /// Position::legal() tests whether a pseudo-legal move is legal

    bool Position::legal(Move m) const {

        assert(is_ok(m));

        Color us = sideToMove;
        Square from = from_sq(m);
        Square to = to_sq(m);

        assert(color_of(moved_piece(m)) == us);
        assert(piece_on(square<KING>(us)) == make_piece(us, KING));
        assert(board_bb() & to);

        Bitboard occupied = (pieces() ^ from) | to;

        // Flying general rule
        Square s = type_of(moved_piece(m)) == KING ? to : square<KING>(us);
        if (attacks_bb(~us, ROOK, s, occupied) & pieces(~us, KING) & ~square_bb(to))
            return false;

        // If the moving piece is a king, check whether the destination square is
        // attacked by the opponent.
        if (type_of(moved_piece(m)) == KING)
            return !attackers_to(to, occupied, ~us);

        // A non-king move is legal if the king is not under attack after the move.
        return !(attackers_to(square<KING>(us), occupied, ~us) & ~SquareBB[to]);
    }


    /// Position::pseudo_legal() takes a random move and tests whether the move is
    /// pseudo legal. It is used to validate moves from TT that can be corrupted
    /// due to SMP concurrent access or hash position key aliasing.

    bool Position::pseudo_legal(const Move m) const {

        Color us = sideToMove;
        Square from = from_sq(m);
        Square to = to_sq(m);
        Piece pc = moved_piece(m);

        // Illegal moves to squares outside of board
        if (!(board_bb() & to))
            return false;

        // If the 'from' square is not occupied by a piece belonging to the side to
        // move, the move is obviously not legal.
        if (pc == NO_PIECE || color_of(pc) != us)
            return false;

        // The destination square cannot be occupied by a friendly piece
        if (pieces(us) & to)
            return false;

        if (!((capture(m) ? attacks_from(us, type_of(pc), from) : moves_from(us, type_of(pc), from)) & to))
            return false;

        // Evasions generator already takes care to avoid some kind of illegal moves
        // and legal() relies on this. We therefore have to take care that the same
        // kind of moves are filtered out here.
        if (checkers() && !(checkers() & non_sliding_riders()))
        {
            if (type_of(pc) != KING)
            {
                // Double check? In this case a king move is required
                if (more_than_one(checkers()))
                    return false;

                // Our move must be a blocking evasion or a capture of the checking piece
                Square checksq = lsb(checkers());
                if (!(between_bb(square<KING>(us), lsb(checkers())) & to)
                    || ((LeaperAttacks[~us][type_of(piece_on(checksq))][checksq] & square<KING>(us)) && !(checkers() & to)))
                    return false;
            }
            // In case of king moves under check we have to remove king so as to catch
            // invalid moves like b1a1 when opposite queen is on c1.
            else if (attackers_to(to, pieces() ^ from, ~us))
                return false;
        }

        return true;
    }


    /// Position::gives_check() tests whether a pseudo-legal move gives a check

    bool Position::gives_check(Move m) const {

        assert(is_ok(m));
        assert(color_of(moved_piece(m)) == sideToMove);

        Square from = from_sq(m);
        Square to = to_sq(m);

        // Is there a direct check?
        PieceType pt = type_of(moved_piece(m));
        if (AttackRiderTypes[pt] & (HOPPING_RIDERS | ASYMMETRICAL_RIDERS))
        {
            Bitboard occupied = (pieces() ^ from) | to;
            if (attacks_bb(sideToMove, pt, to, occupied) & square<KING>(~sideToMove))
                return true;
        }
        else if (check_squares(pt) & to)
            return true;

        // Is there a discovered check?
        if (((blockers_for_king(~sideToMove) & from)) || (non_sliding_riders() & pieces(sideToMove))
            && attackers_to(square<KING>(~sideToMove), (pieces() ^ from) | to, sideToMove))
            return true;

        return false;
    }

    /// Position::do_move() makes a move, and saves all information necessary
    /// to a StateInfo object. The move is assumed to be legal. Pseudo-legal
    /// moves should be filtered out before this function is called.

    void Position::do_move(Move m, StateInfo& newSt, bool givesCheck) {

        assert(is_ok(m));
        assert(&newSt != st);

#ifndef NO_THREADS
        thisThread->nodes.fetch_add(1, std::memory_order_relaxed);
#endif
        Key k = st->key ^ Zobrist::side;

        // Copy some fields of the old state to our new StateInfo object except the
        // ones which are going to be recalculated from scratch anyway and then switch
        // our state pointer to point to the new (ready to be updated) state.
        std::memcpy(static_cast<void*>(&newSt), static_cast<void*>(st), offsetof(StateInfo, key));
        newSt.previous = st;
        st = &newSt;
        st->move = m;

        // Increment ply counters. In particular, rule50 will be reset to zero later on
        // in case of a capture or a pawn move.
        ++gamePly;
        ++st->rule50;
        ++st->pliesFromNull;

        // Used by NNUE
        st->accumulator.computed[WHITE] = false;
        st->accumulator.computed[BLACK] = false;
        auto& dp = st->dirtyPiece;
        dp.dirty_num = 1;

        Color us = sideToMove;
        Color them = ~us;
        Square from = from_sq(m);
        Square to = to_sq(m);
        Piece pc = moved_piece(m);
        Piece captured = piece_on(to);
        st->capturedpromoted = is_promoted(to);
        st->unpromotedCapturedPiece = captured ? unpromoted_piece_on(to) : NO_PIECE;
        st->pass = false;

        assert(color_of(pc) == us);
        assert(captured == NO_PIECE || color_of(captured) == them);
        assert(type_of(captured) != KING);

        if (captured)
        {
            Square capsq = to;

            // update non-pawn material.
            st->nonPawnMaterial[them] -= PieceValue[MG][captured];

            if (Eval::useNNUE)
            {
                dp.dirty_num = 2;  // 1 piece moved, 1 piece captured
                dp.piece[1] = captured;
                dp.from[1] = capsq;
                dp.to[1] = SQ_NONE;
            }

            // Update board and piece lists
            bool capturedPromoted = is_promoted(capsq);
            Piece unpromotedCaptured = unpromoted_piece_on(capsq);
            remove_piece(capsq);

            if (Eval::useNNUE)
                dp.handPiece[1] = NO_PIECE;

            // Update material hash key and prefetch access to materialTable
            k ^= Zobrist::psq[captured][capsq];
            st->materialKey ^= Zobrist::psq[captured][pieceCount[captured]];
#ifndef NO_THREADS
            prefetch(thisThread->materialTable[st->materialKey]);
#endif
            // Reset rule 50 counter
            st->rule50 = 0;
        }
        k ^= Zobrist::psq[pc][from] ^ Zobrist::psq[pc][to];

        // Move the piece.
        if (Eval::useNNUE)
        {
            dp.piece[0] = pc;
            dp.from[0] = from;
            dp.to[0] = to;
        }

        move_piece(from, to);

        // Set capture piece
        st->capturedPiece = captured;

        // Update the key with the final value
        st->key = k;
        // Calculate checkers bitboard (if move gives check)
        st->checkersBB = givesCheck ? attackers_to(square<KING>(them), us) & pieces(us) : Bitboard(0);

        sideToMove = ~sideToMove;

        // Update king attacks used for fast check detection
        set_check_info(st);

        // Calculate the repetition info. It is the ply distance from the previous
        // occurrence of the same position, negative in the 3-fold case, or zero
        // if the position was not repeated.
        st->repetition = 0;
        int end = std::min(st->rule50, st->pliesFromNull);
        if (end >= 4)
        {
            StateInfo* stp = st->previous->previous;
            for (int i = 4; i <= end; i += 2)
            {
                stp = stp->previous->previous;
                if (stp->key == st->key)
                {
                    st->repetition = stp->repetition ? -i : i;
                    break;
                }
            }
        }

        assert(pos_is_ok());
    }


    /// Position::undo_move() unmakes a move. When it returns, the position should
    /// be restored to exactly the same state as before the move was made.

    void Position::undo_move(Move m) {

        assert(is_ok(m));

        sideToMove = ~sideToMove;

        Color us = sideToMove;
        Square from = from_sq(m);
        Square to = to_sq(m);
        Piece pc = piece_on(to);

        assert(empty(from));
        assert(type_of(st->capturedPiece) != KING);

        move_piece(to, from); // Put the piece back at the source square

        if (st->capturedPiece)
        {
            Square capsq = to;

            put_piece(st->capturedPiece, capsq, st->capturedpromoted, st->unpromotedCapturedPiece); // Restore the captured piece
        }

        // Finally point our state pointer back to the previous state
        st = st->previous;
        --gamePly;

        assert(pos_is_ok());
    }


    /// Position::do_null_move() is used to do a "null move": it flips
    /// the side to move without executing any move on the board.

    void Position::do_null_move(StateInfo& newSt) {

        assert(!checkers());
        assert(&newSt != st);

        std::memcpy(&newSt, st, offsetof(StateInfo, accumulator));

        newSt.previous = st;
        st = &newSt;

        st->dirtyPiece.dirty_num = 0;
        st->dirtyPiece.piece[0] = NO_PIECE; // Avoid checks in UpdateAccumulator()
        st->accumulator.computed[WHITE] = false;
        st->accumulator.computed[BLACK] = false;

        st->key ^= Zobrist::side;
        prefetch(TT.first_entry(key()));

        ++st->rule50;
        st->pliesFromNull = 0;

        sideToMove = ~sideToMove;

        set_check_info(st);

        st->repetition = 0;

        assert(pos_is_ok());
    }


    /// Position::undo_null_move() must be used to undo a "null move"

    void Position::undo_null_move() {

        assert(!checkers());

        st = st->previous;
        sideToMove = ~sideToMove;
    }


    /// Position::key_after() computes the new hash key after the given move. Needed
    /// for speculative prefetch. It doesn't recognize special moves like castling,
    /// en passant and promotions.

    Key Position::key_after(Move m) const {

        Square from = from_sq(m);
        Square to = to_sq(m);
        Piece pc = moved_piece(m);
        Piece captured = piece_on(to);
        Key k = st->key ^ Zobrist::side;

        if (captured)
        {
            k ^= Zobrist::psq[captured][to];
        }

        return k ^ Zobrist::psq[pc][to] ^ Zobrist::psq[pc][from];
    }


    Value Position::blast_see(Move m) const {
        assert(is_ok(m));

        Square from = from_sq(m);
        Square to = to_sq(m);
        Color us = color_of(moved_piece(m));
        Bitboard fromto = from | to;
        Bitboard blast = ((attacks_bb<KING>(to) & ~pieces(PAWN)) | fromto) & pieces();

        Value result = VALUE_ZERO;

        // Add the least valuable attacker for quiet moves
        if (!capture(m))
        {
            Bitboard attackers = attackers_to(to, pieces() ^ fromto, ~us);
            Value minAttacker = VALUE_INFINITE;

            if (minAttacker == VALUE_INFINITE)
                return VALUE_ZERO;

            result += minAttacker;
        }

        // Sum up blast piece values
        while (blast)
        {
            Piece bpc = piece_on(pop_lsb(blast));
            result += color_of(bpc) == us ? -CapturePieceValue[MG][bpc] : CapturePieceValue[MG][bpc];
        }

        return capture(m) ? result - 1 : std::min(result, VALUE_ZERO);
    }


    /// Position::see_ge (Static Exchange Evaluation Greater or Equal) tests if the
    /// SEE value of move is greater or equal to the given threshold. We'll use an
    /// algorithm similar to alpha-beta pruning with a null window.

    bool Position::see_ge(Move m, Value threshold) const {

        assert(is_ok(m));

        Square from = from_sq(m), to = to_sq(m);

        int swap = PieceValue[MG][piece_on(to)] - threshold;
        if (swap < 0)
            return false;

        swap = PieceValue[MG][moved_piece(m)] - swap;
        if (swap <= 0)
            return true;

        Bitboard occupied = (pieces() ^ from) ^ to;
        Color stm = color_of(moved_piece(m));
        Bitboard attackers = attackers_to(to, occupied);
        Bitboard stmAttackers, bb;
        int res = 1;
        if (attackers & pieces(stm, KING))
            attackers |= attacks_bb(stm, ROOK, to, occupied & ~pieces(ROOK)) & pieces(~stm, KING);
        if (attackers & pieces(~stm, KING))
            attackers |= attacks_bb(~stm, ROOK, to, occupied & ~pieces(ROOK)) & pieces(stm, KING);


        while (true)
        {
            stm = ~stm;
            attackers &= occupied;

            // If stm has no more attackers then give up: stm loses
            if (!(stmAttackers = attackers & pieces(stm)))
                break;

            // Don't allow pinned pieces to attack as long as there are
            // pinners on their original square.
            if (pinners(~stm) & occupied)
                stmAttackers &= ~blockers_for_king(stm);

            if (!stmAttackers)
                break;

            res ^= 1;

            // Locate and remove the next least valuable attacker, and add to
            // the bitboard 'attackers' any X-ray attackers behind it.
            if ((bb = stmAttackers & pieces(ROOK)))
            {
                if ((swap = RookValueMg - swap) < res)
                    break;

                occupied ^= least_significant_square_bb(bb);
                attackers |= attacks_bb<ROOK>(to, occupied) & pieces(ROOK);
            }

            // fairy pieces
            // pick next piece without considering value
            else if ((bb = stmAttackers & ~pieces(KING)))
            {
                if ((swap = PieceValue[MG][piece_on(lsb(bb))] - swap) < res)
                    break;

                occupied ^= lsb(bb);
            }

            else // KING
                 // If we "capture" with the king but opponent still has attackers,
                 // reverse the result.
                return (attackers & ~pieces(stm)) ? res ^ 1 : res;
        }

        return bool(res);
    }

    /// Position::is_optinal_game_end() tests whether the position may end the game by
    /// 50-move rule, by repetition, or a variant rule that allows a player to claim a game result.

    bool Position::is_optional_game_end(Value& result, int ply, int countStarted) const {

        // n-fold repetition
        if (n_fold_rule())
        {
            int end = std::min(st->rule50, st->pliesFromNull);

            if (end >= 4)
            {
                StateInfo* stp = st->previous->previous;
                int cnt = 0;
                bool perpetualThem = st->checkersBB && stp->checkersBB;
                bool perpetualUs = st->previous->checkersBB && stp->previous->checkersBB;
                Bitboard chaseThem = undo_move_board(st->chased, st->previous->move) & stp->chased;
                Bitboard chaseUs = undo_move_board(st->previous->chased, stp->move) & stp->previous->chased;

                for (int i = 4; i <= end; i += 2)
                {
                    // Chased pieces are empty when there is no previous move
                    if (i != st->pliesFromNull)
                        chaseThem = undo_move_board(chaseThem, stp->previous->move) & stp->previous->previous->chased;
                    stp = stp->previous->previous;
                    perpetualThem &= bool(stp->checkersBB);

                    // Return a draw score if a position repeats once earlier but strictly
                    // after the root, or repeats twice before or at the root.
                    if (stp->key == st->key
                        && ++cnt + 1 == (ply > i ? 2 : n_fold_rule()))
                    {
                        result = convert_mate_value((perpetualThem || perpetualUs) ? (!perpetualUs ? VALUE_MATE : !perpetualThem ? -VALUE_MATE : VALUE_DRAW)
                            : (chaseThem || chaseUs) ? (!chaseUs ? VALUE_MATE : !chaseThem ? -VALUE_MATE : VALUE_DRAW)
                            : VALUE_DRAW, ply);
                        return true;
                    }

                    if (i + 1 <= end)
                    {
                        perpetualUs &= bool(stp->previous->checkersBB);
                        chaseUs = undo_move_board(chaseUs, stp->move) & stp->previous->chased;
                    }
                }
            }
        }

        return false;
    }

    /// Position::is_immediate_game_end() tests whether the position ends the game
    /// immediately by a variant rule, i.e., there are no more legal moves.
    /// It does not not detect stalemates.

    bool Position::is_immediate_game_end(Value& result, int ply) const {
        return false;
    }

    // Position::chased() tests whether the last move was a chase.

    Bitboard Position::chased() const {
        Bitboard b = 0;
        if (st->move == MOVE_NONE)
            return b;

        Bitboard pins = blockers_for_king(sideToMove);
        Square ourKing = square<KING>(sideToMove);
        Square oppKing = square<KING>(~sideToMove);
        if (file_bb(file_of(ourKing)) & file_bb(file_of(oppKing))) {
            Bitboard kingFilePieces = between_bb(ourKing, oppKing) ^ square_bb(oppKing);
            if (!more_than_one(kingFilePieces & pieces()))
                pins |= kingFilePieces & pieces(sideToMove);
        }
        auto addChased = [&](Square attackerSq, PieceType attackerType, Bitboard attacks) {
            if (attacks & ~b)
            {
                // Exclude attacks on unpromoted soldiers and checks
                attacks &= ~(pieces(sideToMove, KING, SOLDIER) ^ promoted_soldiers(sideToMove));
                // Attacks against stronger pieces
                if (attackerType == HORSE || attackerType == CANNON)
                    b |= attacks & pieces(sideToMove, ROOK);
                if (attackerType == ELEPHANT || attackerType == FERS)
                    b |= attacks & pieces(sideToMove, ROOK, CANNON, HORSE);
                // Exclude mutual/symmetric attacks
                // Exceptions:
                // - asymmetric pieces ("impaired horse")
                // - pins
                if (attackerType == HORSE && (PseudoAttacks[WHITE][FERS][attackerSq] & pieces()))
                {
                    Bitboard horses = attacks & pieces(sideToMove, attackerType);
                    while (horses)
                    {
                        Square s = pop_lsb(horses);
                        if (attacks_bb(sideToMove, attackerType, s, pieces()) & attackerSq)
                            attacks ^= s;
                    }
                }
                else
                    attacks &= ~pieces(sideToMove, attackerType) | pins;
                // Attacks against potentially unprotected pieces
                while (attacks)
                {
                    Square s = pop_lsb(attacks);
                    Bitboard roots = attackers_to(s, pieces() ^ attackerSq, sideToMove) & ~pins;
                    if (!roots || (roots == pieces(sideToMove, KING) && (attacks_bb(sideToMove, ROOK, square<KING>(~sideToMove), pieces() ^ attackerSq) & s)))
                        b |= s;
                }
            }
        };

        // Direct attacks
        Square from = from_sq(st->move);
        Square to = to_sq(st->move);
        PieceType movedPiece = type_of(piece_on(to));
        if (movedPiece != KING && movedPiece != SOLDIER)
        {
            Bitboard directAttacks = attacks_from(~sideToMove, movedPiece, to) & pieces(sideToMove);
            // Only new attacks count. This avoids expensive comparison of previous and new attacks.
            if (movedPiece == ROOK || movedPiece == CANNON)
                directAttacks &= ~line_bb(from, to);
            addChased(to, movedPiece, directAttacks);
        }

        // Discovered attacks
        Bitboard discoveryCandidates = (PseudoAttacks[WHITE][WAZIR][from] & pieces(~sideToMove, HORSE))
            | (PseudoAttacks[WHITE][FERS][from] & pieces(~sideToMove, ELEPHANT))
            | (PseudoAttacks[WHITE][ROOK][from] & pieces(~sideToMove, CANNON, ROOK))
            | (PseudoAttacks[WHITE][ROOK][to] & pieces(~sideToMove, CANNON));
        while (discoveryCandidates)
        {
            Square s = pop_lsb(discoveryCandidates);
            PieceType discoveryPiece = type_of(piece_on(s));
            Bitboard discoveries = pieces(sideToMove)
                & attacks_bb(~sideToMove, discoveryPiece, s, pieces())
                & ~attacks_bb(~sideToMove, discoveryPiece, s, (captured_piece() ? pieces() : pieces() ^ to) ^ from);
            addChased(s, discoveryPiece, discoveries);
        }

        // Changes in real roots and discovered checks
        if (st->pliesFromNull > 0)
        {
            // Fake roots
            Bitboard newPins = st->blockersForKing[sideToMove] & ~st->previous->blockersForKing[sideToMove] & pieces(sideToMove);
            while (newPins)
            {
                Square s = pop_lsb(newPins);
                PieceType pinnedPiece = type_of(piece_on(s));
                Bitboard fakeRooted = pieces(sideToMove)
                    & ~(pieces(sideToMove, KING, SOLDIER) ^ promoted_soldiers(sideToMove))
                    & attacks_bb(sideToMove, pinnedPiece, s, pieces());
                while (fakeRooted)
                {
                    Square s2 = pop_lsb(fakeRooted);
                    if (attackers_to(s2, ~sideToMove) & ~blockers_for_king(~sideToMove))
                        b |= s2;
                }
            }
            // Discovered checks
            Bitboard newDiscoverers = st->blockersForKing[sideToMove] & ~st->previous->blockersForKing[sideToMove] & pieces(~sideToMove);
            while (newDiscoverers)
            {
                Square s = pop_lsb(newDiscoverers);
                PieceType discoveryPiece = type_of(piece_on(s));
                Bitboard discoveryAttacks = attacks_from(~sideToMove, discoveryPiece, s) & pieces(sideToMove);
                // Include all captures except where the king can pseudo-legally recapture
                b |= discoveryAttacks & ~attacks_from(sideToMove, KING, square<KING>(sideToMove));
                // Include captures where king can not legally recapture
                discoveryAttacks &= attacks_from(sideToMove, KING, square<KING>(sideToMove));
                while (discoveryAttacks)
                {
                    Square s2 = pop_lsb(discoveryAttacks);
                    if (attackers_to(s2, pieces() ^ s ^ square<KING>(sideToMove), ~sideToMove) & ~square_bb(s))
                        b |= s2;
                }
            }
        }

        return b;
    }

    // Position::has_repeated() tests whether there has been at least one repetition
    // of positions since the last capture or pawn move.

    bool Position::has_repeated() const {

        StateInfo* stc = st;
        int end = std::min(st->rule50, st->pliesFromNull);
        while (end-- >= 4)
        {
            if (stc->repetition)
                return true;

            stc = stc->previous;
        }
        return false;
    }


    /// Position::flip() flips position with the white and black sides reversed. This
    /// is only useful for debugging e.g. for finding evaluation symmetry bugs.

    void Position::flip() {

        string f, token;
        std::stringstream ss(fen());

        for (Rank r = max_rank(); r >= RANK_1; --r) // Piece placement
        {
            std::getline(ss, token, r > RANK_1 ? '/' : ' ');
            f.insert(0, token + (f.empty() ? " " : "/"));
        }

        ss >> token; // Active color
        f += (token == "w" ? "B " : "W "); // Will be lowercased later

        ss >> token; // Castling availability
        f += token + " ";

        std::transform(f.begin(), f.end(), f.begin(),
            [](char c) { return char(islower(c) ? toupper(c) : tolower(c)); });

        ss >> token; // En passant square
        f += (token == "-" ? token : token.replace(1, 1, token[1] == '3' ? "6" : "3"));

        std::getline(ss, token); // Half and full moves
        f += token;

        set(variant(), f, st, this_thread());

        assert(pos_is_ok());
    }


    /// Position::pos_is_ok() performs some consistency checks for the
    /// position object and raises an asserts if something wrong is detected.
    /// This is meant to be helpful when debugging.

    bool Position::pos_is_ok() const {

        constexpr bool Fast = true; // Quick (default) or full check?

        if ((sideToMove != WHITE && sideToMove != BLACK)
            || (piece_on(square<KING>(WHITE)) != make_piece(WHITE, KING))
            || (piece_on(square<KING>(BLACK)) != make_piece(BLACK, KING)))
            assert(0 && "pos_is_ok: Default");

        if (Fast)
            return true;

        if (pieceCount[make_piece(~sideToMove, KING)]
            && (attackers_to(square<KING>(~sideToMove)) & pieces(sideToMove)))
            assert(0 && "pos_is_ok: Kings");

        if (pieceCount[make_piece(WHITE, PAWN)] > 64
            || pieceCount[make_piece(BLACK, PAWN)] > 64)
            assert(0 && "pos_is_ok: Pawns");

        if ((pieces(WHITE) & pieces(BLACK))
            || (pieces(WHITE) | pieces(BLACK)) != pieces()
            || popcount(pieces(WHITE)) > 64
            || popcount(pieces(BLACK)) > 64)
            assert(0 && "pos_is_ok: Bitboards");

        for (PieceType p1 = PAWN; p1 <= KING; ++p1)
            for (PieceType p2 = PAWN; p2 <= KING; ++p2)
                if (p1 != p2 && (pieces(p1) & pieces(p2)))
                    assert(0 && "pos_is_ok: Bitboards");

        StateInfo si = *st;
        ASSERT_ALIGNED(&si, Eval::NNUE::CacheLineSize);

        set_state(&si);
        if (std::memcmp(&si, st, sizeof(StateInfo)))
            assert(0 && "pos_is_ok: State");

        for (Color c : {WHITE, BLACK})
            for (PieceType pt = PAWN; pt <= KING; ++pt)
            {
                Piece pc = make_piece(c, pt);
                if (pieceCount[pc] != popcount(pieces(c, pt))
                    || pieceCount[pc] != std::count(board, board + SQUARE_NB, pc))
                    assert(0 && "pos_is_ok: Pieces");
            }

        return true;
    }

} // namespace Stockfish

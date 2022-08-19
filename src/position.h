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

#ifndef POSITION_H_INCLUDED
#define POSITION_H_INCLUDED

#include <cassert>
#include <deque>
#include <memory> // For std::unique_ptr
#include <string>
#include <functional>

#include "bitboard.h"
#include "evaluate.h"
#include "psqt.h"
#include "types.h"
#include "variant.h"
#include "movegen.h"

#include "nnue/nnue_accumulator.h"

namespace Stockfish {

/// StateInfo struct stores information needed to restore a Position object to
/// its previous state when we retract a move. Whenever a move is made on the
/// board (by calling Position::do_move), a StateInfo object must be passed.

struct StateInfo {

  // Copied when making a move
  Key    pawnKey;
  Key    materialKey;
  Value  nonPawnMaterial[COLOR_NB];
  int    rule50;
  int    pliesFromNull;
  int    countingPly;
  CheckCount checksRemaining[COLOR_NB];
  Square epSquare;
  Bitboard gatesBB[COLOR_NB];

  // Not copied when making a move (will be recomputed anyhow)
  Key        key;
  Bitboard   checkersBB;
  Piece      unpromotedCapturedPiece;
  Piece      unpromotedBycatch[SQUARE_NB];
  Bitboard   promotedBycatch;
  Bitboard   demotedBycatch;
  StateInfo* previous;
  Bitboard   blockersForKing[COLOR_NB];
  Bitboard   pinners[COLOR_NB];
  Bitboard   checkSquares[PIECE_TYPE_NB];
  Piece      capturedPiece;
  Bitboard   nonSlidingRiders;
  Bitboard   pseudoRoyals;
  OptBool    legalCapture;
  bool       capturedpromoted;
  bool       shak;
  Bitboard   chased;
  bool       pass;
  Move       move;
  int        repetition;

  // Used by NNUE
  Eval::NNUE::Accumulator accumulator;
  DirtyPiece dirtyPiece;
};


/// A list to keep track of the position states along the setup moves (from the
/// start position to the position just before the search starts). Needed by
/// 'draw by repetition' detection. Use a std::deque because pointers to
/// elements are not invalidated upon list resizing.
typedef std::unique_ptr<std::deque<StateInfo>> StateListPtr;


/// Position class stores information regarding the board representation as
/// pieces, side to move, hash keys, castling info, etc. Important methods are
/// do_move() and undo_move(), used by the search to update node info when
/// traversing the search tree.
class Thread;

class Position {
public:
  static void init();

  Position() = default;
  Position(const Position&) = delete;
  Position& operator=(const Position&) = delete;

  // FEN string input/output
  Position& set(const Variant* v, const std::string& fenStr, StateInfo* si, Thread* th, bool sfen = false);
  Position& set(const std::string& code, Color c, StateInfo* si);
  std::string fen(bool sfen = false, bool showPromoted = false, int countStarted = 0, std::string holdings = "-") const;

  // Variant rule properties
  const Variant* variant() const;
  Rank max_rank() const;
  File max_file() const;
  int ranks() const;
  int files() const;
  Bitboard board_bb() const;
  Bitboard board_bb(Color c, PieceType pt) const;
  const std::string& piece_to_char() const;
  const std::string& piece_to_char_synonyms() const;
  Rank promotion_rank() const;
  PieceType nnue_king() const;
  Square nnue_king_square(Color c) const;
  bool nnue_applicable() const;
  bool has_capture() const;
  bool can_drop(Color c, PieceType pt) const;
  Bitboard drop_region(Color c, PieceType pt) const;
  Bitboard promoted_soldiers(Color c) const;
  // winning conditions
  int n_fold_rule() const;
  Value stalemate_value(int ply = 0) const;
  Value checkmate_value(int ply = 0) const;
  CheckCount checks_remaining(Color c) const;

  // Variant-specific properties
  int count_in_hand(PieceType pt) const;
  int count_in_hand(Color c, PieceType pt) const;
  int count_with_hand(Color c, PieceType pt) const;

  // Position representation
  Bitboard pieces(PieceType pt = ALL_PIECES) const;
  Bitboard pieces(PieceType pt1, PieceType pt2) const;
  Bitboard pieces(Color c) const;
  Bitboard pieces(Color c, PieceType pt) const;
  Bitboard pieces(Color c, PieceType pt1, PieceType pt2) const;
  Bitboard pieces(Color c, PieceType pt1, PieceType pt2, PieceType pt3) const;
  Bitboard non_sliding_riders() const;
  Piece piece_on(Square s) const;
  Piece unpromoted_piece_on(Square s) const;
  Square ep_square() const;
  Bitboard gates(Color c) const;
  bool empty(Square s) const;
  int count(Color c, PieceType pt) const;
  template<PieceType Pt> int count(Color c) const;
  template<PieceType Pt> int count() const;
  template<PieceType Pt> Square square(Color c) const;
  Square square(Color c, PieceType pt) const;
  bool is_on_semiopen_file(Color c, Square s) const;

  // Checking
  Bitboard checkers() const;
  Bitboard blockers_for_king(Color c) const;
  Bitboard check_squares(PieceType pt) const;
  Bitboard pinners(Color c) const;
  Bitboard attackers_to_pseudo_royals(Color c) const;

  // Attacks to/from a given square
  Bitboard attackers_to(Square s) const;
  Bitboard attackers_to(Square s, Color c) const;
  Bitboard attackers_to(Square s, Bitboard occupied) const;
  Bitboard attackers_to(Square s, Bitboard occupied, Color c) const;
  Bitboard attacks_from(Color c, PieceType pt, Square s) const;
  Bitboard moves_from(Color c, PieceType pt, Square s) const;
  Bitboard slider_blockers(Bitboard sliders, Square s, Bitboard& pinners, Color c) const;

  // Properties of moves
  bool legal(Move m) const;
  bool pseudo_legal(const Move m) const;
  bool capture(Move m) const;
  bool capture_or_promotion(Move m) const;
  bool gives_check(Move m) const;
  Piece moved_piece(Move m) const;
  Piece captured_piece() const;

  // Piece specific
  bool is_promoted(Square s) const;

  // Doing and undoing moves
  void do_move(Move m, StateInfo& newSt);
  void do_move(Move m, StateInfo& newSt, bool givesCheck);
  void undo_move(Move m);
  void do_null_move(StateInfo& newSt);
  void undo_null_move();

  // Static Exchange Evaluation
  Value blast_see(Move m) const;
  bool see_ge(Move m, Value threshold = VALUE_ZERO) const;

  // Accessing hash keys
  Key key() const;
  Key key_after(Move m) const;
  Key material_key() const;
  Key pawn_key() const;

  // Other properties of the position
  Color side_to_move() const;
  int game_ply() const;
  Thread* this_thread() const;
  bool is_immediate_game_end() const;
  bool is_immediate_game_end(Value& result, int ply = 0) const;
  bool is_optional_game_end() const;
  bool is_optional_game_end(Value& result, int ply = 0, int countStarted = 0) const;
  bool is_game_end(Value& result, int ply = 0) const;
  bool is_draw(int ply) const;
  bool has_repeated() const;
  Bitboard chased() const;
  int counting_ply(int countStarted) const;
  int rule50_count() const;
  Score psq_score() const;
  Value psq_eg_stm() const;
  Value non_pawn_material(Color c) const;
  Value non_pawn_material() const;

  // Position consistency check, for debugging
  bool pos_is_ok() const;
  void flip();

  // Used by NNUE
  StateInfo* state() const;

  void put_piece(Piece pc, Square s, bool isPromoted = false, Piece unpromotedPc = NO_PIECE);
  void remove_piece(Square s);

private:
  // Initialization helpers (used while setting up a position)
  void set_state(StateInfo* si) const;
  void set_check_info(StateInfo* si) const;

  // Other helpers
  void move_piece(Square from, Square to);

  // Data members
  Piece board[SQUARE_NB];
  Piece unpromotedBoard[SQUARE_NB];
  Bitboard byTypeBB[PIECE_TYPE_NB];
  Bitboard byColorBB[COLOR_NB];
  int pieceCount[PIECE_NB];
  Thread* thisThread;
  StateInfo* st;
  int gamePly;
  Color sideToMove;
  Score psq;

  // variant-specific
  const Variant* var;
  int pieceCountInHand[COLOR_NB][PIECE_TYPE_NB];
  Bitboard promotedPieces;
  void add_to_hand(Piece pc);
  void remove_from_hand(Piece pc);
};

extern std::ostream& operator<<(std::ostream& os, const Position& pos);

inline const Variant* Position::variant() const {
  assert(var != nullptr);
  return var;
}

inline Rank Position::max_rank() const {
  assert(var != nullptr);
  return var->maxRank;
}

inline File Position::max_file() const {
  assert(var != nullptr);
  return var->maxFile;
}

inline int Position::ranks() const {
  assert(var != nullptr);
  return var->maxRank + 1;
}

inline int Position::files() const {
  assert(var != nullptr);
  return var->maxFile + 1;
}

inline Bitboard Position::board_bb() const {
  assert(var != nullptr);
  return board_size_bb(var->maxFile, var->maxRank);
}

inline Bitboard Position::board_bb(Color c, PieceType pt) const {
  return BoardBB[c][pt];
}

inline const std::string& Position::piece_to_char() const {
  assert(var != nullptr);
  return var->pieceToChar;
}

inline const std::string& Position::piece_to_char_synonyms() const {
  assert(var != nullptr);
  return var->pieceToCharSynonyms;
}

inline Rank Position::promotion_rank() const {
  assert(var != nullptr);
  return RANK_8;
}

inline PieceType Position::nnue_king() const {
  assert(var != nullptr);
  return KING;
}

inline Square Position::nnue_king_square(Color c) const {
  return nnue_king() ? square(c, nnue_king()) : SQ_NONE;
}

inline bool Position::nnue_applicable() const {
  // Do not use NNUE during setup phases (placement, sittuyin)
  return !count_in_hand(ALL_PIECES);
}

inline bool Position::has_capture() const {
  // Check for cached value
  if (st->legalCapture != NO_VALUE)
      return st->legalCapture == VALUE_TRUE;
  if (checkers())
  {
      for (const auto& mevasion : MoveList<EVASIONS>(*this))
          if (capture(mevasion) && legal(mevasion))
          {
              st->legalCapture = VALUE_TRUE;
              return true;
          }
  }
  else
  {
      for (const auto& mcap : MoveList<CAPTURES>(*this))
          if (capture(mcap) && legal(mcap))
          {
              st->legalCapture = VALUE_TRUE;
              return true;
          }
  }
  st->legalCapture = VALUE_FALSE;
  return false;
}


inline Bitboard Position::drop_region(Color c, PieceType pt) const {
  Bitboard b = AllSquares & board_bb(c, pt);

  // Pawns on back ranks
  if (pt == PAWN)
  {
      b &= ~zone_bb(c, promotion_rank(), max_rank());
      b &= ~rank_bb(relative_rank(c, RANK_1, max_rank()));
  }

  return b;
}

inline Bitboard Position::promoted_soldiers(Color c) const {
  assert(var != nullptr);
  return pieces(c, SOLDIER) & zone_bb(c, RANK_6, max_rank());
}

inline int Position::n_fold_rule() const {
  assert(var != nullptr);
  return 3;
}

inline Value Position::stalemate_value(int ply) const {
  assert(var != nullptr);
  return convert_mate_value(-VALUE_MATE, ply);
}

inline Value Position::checkmate_value(int ply) const {
  assert(var != nullptr);
  // Return mate value
  return convert_mate_value(-VALUE_MATE, ply);
}

inline CheckCount Position::checks_remaining(Color c) const {
  return st->checksRemaining[c];
}

inline bool Position::is_immediate_game_end() const {
  Value result;
  return is_immediate_game_end(result);
}

inline bool Position::is_optional_game_end() const {
  Value result;
  return is_optional_game_end(result);
}

inline bool Position::is_draw(int ply) const {
  Value result;
  return is_optional_game_end(result, ply);
}

inline bool Position::is_game_end(Value& result, int ply) const {
  return is_immediate_game_end(result, ply) || is_optional_game_end(result, ply);
}

inline Color Position::side_to_move() const {
  return sideToMove;
}

inline Piece Position::piece_on(Square s) const {
  assert(is_ok(s));
  return board[s];
}

inline bool Position::empty(Square s) const {
  return piece_on(s) == NO_PIECE;
}

inline Piece Position::unpromoted_piece_on(Square s) const {
  return unpromotedBoard[s];
}

inline Piece Position::moved_piece(Move m) const {
  return piece_on(from_sq(m));
}

inline Bitboard Position::pieces(PieceType pt) const {
  return byTypeBB[pt];
}

inline Bitboard Position::pieces(PieceType pt1, PieceType pt2) const {
  return pieces(pt1) | pieces(pt2);
}

inline Bitboard Position::pieces(Color c) const {
  return byColorBB[c];
}

inline Bitboard Position::pieces(Color c, PieceType pt) const {
  return pieces(c) & pieces(pt);
}

inline Bitboard Position::pieces(Color c, PieceType pt1, PieceType pt2) const {
  return pieces(c) & (pieces(pt1) | pieces(pt2));
}

inline Bitboard Position::pieces(Color c, PieceType pt1, PieceType pt2, PieceType pt3) const {
  return pieces(c) & (pieces(pt1) | pieces(pt2) | pieces(pt3));
}

inline Bitboard Position::non_sliding_riders() const {
  return st->nonSlidingRiders;
}

inline int Position::count(Color c, PieceType pt) const {
  return pieceCount[make_piece(c, pt)];
}

template<PieceType Pt> inline int Position::count(Color c) const {
  return pieceCount[make_piece(c, Pt)];
}

template<PieceType Pt> inline int Position::count() const {
  return count<Pt>(WHITE) + count<Pt>(BLACK);
}

template<PieceType Pt> inline Square Position::square(Color c) const {
  assert(count<Pt>(c) == 1);
  return lsb(pieces(c, Pt));
}

inline Square Position::square(Color c, PieceType pt) const {
  assert(count(c, pt) == 1);
  return lsb(pieces(c, pt));
}

inline Square Position::ep_square() const {
  return st->epSquare;
}

inline Bitboard Position::gates(Color c) const {
  assert(var != nullptr);
  return st->gatesBB[c];
}

inline bool Position::is_on_semiopen_file(Color c, Square s) const {
  return !(pieces(c, SOLDIER) & file_bb(s));
}

inline Bitboard Position::attacks_from(Color c, PieceType pt, Square s) const {

  PieceType movePt = pt == KING ? WAZIR : pt;
  Bitboard b = attacks_bb(c, movePt, s, byTypeBB[ALL_PIECES]);
  // Xiangqi soldier
  if (pt == SOLDIER && !(promoted_soldiers(c) & s))
      b &= file_bb(file_of(s));
  return b & board_bb(c, pt);
}

inline Bitboard Position::moves_from(Color c, PieceType pt, Square s) const {

  PieceType movePt = pt == KING ? WAZIR : pt;
  Bitboard b = moves_bb(c, movePt, s, byTypeBB[ALL_PIECES]);
  // Xiangqi soldier
  if (pt == SOLDIER && !(promoted_soldiers(c) & s))
      b &= file_bb(file_of(s));
  return b & board_bb(c, pt);
}

inline Bitboard Position::attackers_to(Square s) const {
  return attackers_to(s, pieces());
}

inline Bitboard Position::attackers_to(Square s, Color c) const {
  return attackers_to(s, byTypeBB[ALL_PIECES], c);
}

inline Bitboard Position::checkers() const {
  return st->checkersBB;
}

inline Bitboard Position::blockers_for_king(Color c) const {
  return st->blockersForKing[c];
}

inline Bitboard Position::pinners(Color c) const {
  return st->pinners[c];
}

inline Bitboard Position::check_squares(PieceType pt) const {
  return st->checkSquares[pt];
}

inline Key Position::key() const {
  return st->rule50 < 14 ? st->key
                         : st->key ^ make_key((st->rule50 - 14) / 8);
}

inline Key Position::pawn_key() const {
  return st->pawnKey;
}

inline Key Position::material_key() const {
  return st->materialKey;
}

inline Score Position::psq_score() const {
  return psq;
}

inline Value Position::psq_eg_stm() const {
  return (sideToMove == WHITE ? 1 : -1) * eg_value(psq);
}

inline Value Position::non_pawn_material(Color c) const {
  return st->nonPawnMaterial[c];
}

inline Value Position::non_pawn_material() const {
  return non_pawn_material(WHITE) + non_pawn_material(BLACK);
}

inline int Position::game_ply() const {
  return gamePly;
}

inline int Position::counting_ply(int countStarted) const {
  return countStarted == 0 || (count<ALL_PIECES>(WHITE) <= 1 || count<ALL_PIECES>(BLACK) <= 1) ? st->countingPly : countStarted < 0 ? 0 : std::min(st->countingPly, std::max(1 + gamePly - countStarted, 0));
}

inline int Position::rule50_count() const {
  return st->rule50;
}

inline bool Position::is_promoted(Square s) const {
  return promotedPieces & s;
}

inline bool Position::capture_or_promotion(Move m) const {
  assert(is_ok(m));
  return !empty(to_sq(m));
}

inline bool Position::capture(Move m) const {
  assert(is_ok(m));
  // Castling is encoded as "king captures rook"
  return !empty(to_sq(m)) && from_sq(m) != to_sq(m);
}

inline Piece Position::captured_piece() const {
  return st->capturedPiece;
}

inline Thread* Position::this_thread() const {
  return thisThread;
}

inline void Position::put_piece(Piece pc, Square s, bool isPromoted, Piece unpromotedPc) {

  board[s] = pc;
  byTypeBB[ALL_PIECES] |= byTypeBB[type_of(pc)] |= s;
  byColorBB[color_of(pc)] |= s;
  pieceCount[pc]++;
  pieceCount[make_piece(color_of(pc), ALL_PIECES)]++;
  psq += PSQT::psq[pc][s];
  if (isPromoted)
      promotedPieces |= s;
  unpromotedBoard[s] = unpromotedPc;
}

inline void Position::remove_piece(Square s) {

  Piece pc = board[s];
  byTypeBB[ALL_PIECES] ^= s;
  byTypeBB[type_of(pc)] ^= s;
  byColorBB[color_of(pc)] ^= s;
  board[s] = NO_PIECE;
  pieceCount[pc]--;
  pieceCount[make_piece(color_of(pc), ALL_PIECES)]--;
  psq -= PSQT::psq[pc][s];
  promotedPieces -= s;
  unpromotedBoard[s] = NO_PIECE;
}

inline void Position::move_piece(Square from, Square to) {

  Piece pc = board[from];
  Bitboard fromTo = square_bb(from) ^ to; // from == to needs to cancel out
  byTypeBB[ALL_PIECES] ^= fromTo;
  byTypeBB[type_of(pc)] ^= fromTo;
  byColorBB[color_of(pc)] ^= fromTo;
  board[from] = NO_PIECE;
  board[to] = pc;
  psq += PSQT::psq[pc][to] - PSQT::psq[pc][from];
  if (is_promoted(from))
      promotedPieces ^= fromTo;
  unpromotedBoard[to] = unpromotedBoard[from];
  unpromotedBoard[from] = NO_PIECE;
}

inline void Position::do_move(Move m, StateInfo& newSt) {
  do_move(m, newSt, gives_check(m));
}

inline StateInfo* Position::state() const {

  return st;
}

// Variant-specific

inline int Position::count_in_hand(PieceType pt) const {
  return pieceCountInHand[WHITE][pt] + pieceCountInHand[BLACK][pt];
}

inline int Position::count_in_hand(Color c, PieceType pt) const {
  return pieceCountInHand[c][pt];
}

inline int Position::count_with_hand(Color c, PieceType pt) const {
  return pieceCount[make_piece(c, pt)] + pieceCountInHand[c][pt];
}

inline void Position::add_to_hand(Piece pc) {
  pieceCountInHand[color_of(pc)][type_of(pc)]++;
  pieceCountInHand[color_of(pc)][ALL_PIECES]++;
  psq += PSQT::psq[pc][SQ_NONE];
}

inline void Position::remove_from_hand(Piece pc) {
  pieceCountInHand[color_of(pc)][type_of(pc)]--;
  pieceCountInHand[color_of(pc)][ALL_PIECES]--;
  psq -= PSQT::psq[pc][SQ_NONE];
}

inline bool Position::can_drop(Color c, PieceType pt) const {
  return count_in_hand(c, pt) > 0;
}

} // namespace Stockfish

#endif // #ifndef POSITION_H_INCLUDED

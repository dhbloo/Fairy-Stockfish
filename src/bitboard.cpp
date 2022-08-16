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
#include <bitset>

#include "bitboard.h"
#include "magic.h"
#include "misc.h"
#include "piece.h"

namespace Stockfish {

uint8_t PopCnt16[1 << 16];
uint8_t SquareDistance[SQUARE_NB][SQUARE_NB];

Bitboard SquareBB[SQUARE_NB];
Bitboard LineBB[SQUARE_NB][SQUARE_NB];
Bitboard BetweenBB[SQUARE_NB][SQUARE_NB];
Bitboard PseudoAttacks[COLOR_NB][PIECE_TYPE_NB][SQUARE_NB];
Bitboard PseudoMoves[COLOR_NB][PIECE_TYPE_NB][SQUARE_NB];
Bitboard LeaperAttacks[COLOR_NB][PIECE_TYPE_NB][SQUARE_NB];
Bitboard LeaperMoves[COLOR_NB][PIECE_TYPE_NB][SQUARE_NB];
Bitboard BoardSizeBB[FILE_NB][RANK_NB];
RiderType AttackRiderTypes[PIECE_TYPE_NB];
RiderType MoveRiderTypes[PIECE_TYPE_NB];

Magic RookMagicsH[SQUARE_NB];
Magic RookMagicsV[SQUARE_NB];
Magic CannonMagicsH[SQUARE_NB];
Magic CannonMagicsV[SQUARE_NB];
Magic HorseMagics[SQUARE_NB];
Magic ByHorseMagics[SQUARE_NB];
Magic ElephantMagics[SQUARE_NB];

Magic* magics[] = {RookMagicsH, RookMagicsV, CannonMagicsH, CannonMagicsV, HorseMagics, ElephantMagics};

namespace {

// Some magics need to be split in order to reduce memory consumption.
// Otherwise on a 12x10 board they can be >100 MB.
#ifdef LARGEBOARDS
  Bitboard RookTableH[0x11800];  // To store horizontalrook attacks
  Bitboard RookTableV[0x4800];  // To store vertical rook attacks
  Bitboard BishopTable[0x33C00]; // To store bishop attacks
  Bitboard CannonTableH[0x11800];  // To store horizontal cannon attacks
  Bitboard CannonTableV[0x4800];  // To store vertical cannon attacks
  Bitboard HorseTable[0x500];  // To store horse attacks
  Bitboard ByHorseTable[0x1000]; // To store attacks by horse
  Bitboard ElephantTable[0x400];  // To store elephant attacks
#else
  Bitboard RookTableH[0xA00];  // To store horizontal rook attacks
  Bitboard RookTableV[0xA00];  // To store vertical rook attacks
  Bitboard BishopTable[0x1480]; // To store bishop attacks
  Bitboard CannonTableH[0xA00];  // To store horizontal cannon attacks
  Bitboard CannonTableV[0xA00];  // To store vertical cannon attacks
  Bitboard HorseTable[0x240];  // To store horse attacks
  Bitboard ElephantTable[0x1A0];  // To store elephant attacks
#endif

  // Rider directions
  const std::map<Direction, int> RookDirectionsV { {NORTH, 0}, {SOUTH, 0}};
  const std::map<Direction, int> RookDirectionsH { {EAST, 0}, {WEST, 0} };
  const std::map<Direction, int> HorseDirections { {2 * SOUTH + WEST, 0}, {2 * SOUTH + EAST, 0}, {SOUTH + 2 * WEST, 0}, {SOUTH + 2 * EAST, 0},
                                                   {NORTH + 2 * WEST, 0}, {NORTH + 2 * EAST, 0}, {2 * NORTH + WEST, 0}, {2 * NORTH + EAST, 0} };
  const std::map<Direction, int> ElephantDirections { {2 * NORTH_EAST, 0}, {2 * SOUTH_EAST, 0}, {2 * SOUTH_WEST, 0}, {2 * NORTH_WEST, 0} };

  enum MovementType { RIDER, HOPPER, LAME_LEAPER, UNLIMITED_RIDER, BY_HORSE };

  template <MovementType MT>
  void init_magics(Bitboard table[], Magic magics[], std::map<Direction, int> directions, const Bitboard magicsInit[]);

  template <MovementType MT>
  void init_magics(Bitboard table[], Magic magics[], std::map<Direction, int> directions);

  template <MovementType MT>
  Bitboard sliding_attack(std::map<Direction, int> directions, Square sq, Bitboard occupied, Color c = WHITE) {
    assert(MT != LAME_LEAPER);

    Bitboard attack = 0;

    for (auto const& [d, limit] : directions)
    {
        int count = 0;
        bool hurdle = false;
        for (Square s = sq + (c == WHITE ? d : -d);
             is_ok(s) && distance(s, s - (c == WHITE ? d : -d)) <= 2;
             s += (c == WHITE ? d : -d))
        {
            if (MT != HOPPER || hurdle)
            {
                attack |= s;
                if (limit && MT != UNLIMITED_RIDER && ++count >= limit)
                    break;
            }

            if (occupied & s)
            {
                if (MT == HOPPER && !hurdle)
                    hurdle = true;
                else
                    break;
            }
        }
    }

    return attack;
  }

  Bitboard lame_leaper_path(Direction d, Square s) {
    Direction dr = d > 0 ? NORTH : SOUTH;
    Direction df = (std::abs(d % NORTH) < NORTH / 2 ? d % NORTH : -(d % NORTH)) < 0 ? WEST : EAST;
    Square to = s + d;
    Bitboard b = 0;
    if (!is_ok(to) || distance(s, to) >= 4)
        return b;
    while (s != to)
    {
        int diff = std::abs(file_of(to) - file_of(s)) - std::abs(rank_of(to) - rank_of(s));
        if (diff > 0)
            s += df;
        else if (diff < 0)
            s += dr;
        else
            s += df + dr;

        if (s != to)
            b |= s;
    }
    return b;
  }

  Bitboard horse_path(Direction d, Square s) {
    Square to = s + d;
    Bitboard b = 0;
    if (!is_ok(to) || distance(s, to) >= 4)
        return b;

    // use lame_leaper_path but in a reversed order
    return lame_leaper_path(-d, to);
  }

  Bitboard lame_leaper_path(std::map<Direction, int> directions, Square s) {
    Bitboard b = 0;
    for (const auto& i : directions)
        b |= lame_leaper_path(i.first, s);
    return b;
  }

  Bitboard horse_path(std::map<Direction, int> directions, Square s) {
    Bitboard b = 0;
    for (const auto& i : directions)
        b |= horse_path(i.first, s);
    return b;
  }

  Bitboard lame_leaper_attack(std::map<Direction, int> directions, Square s, Bitboard occupied) {
    Bitboard b = 0;
    for (const auto& i : directions)
    {
        Square to = s + i.first;
        if (is_ok(to) && distance(s, to) < 4 && !(lame_leaper_path(i.first, s) & occupied))
            b |= to;
    }
    return b;
  }

  Bitboard horse_attack(std::map<Direction, int> directions, Square s, Bitboard occupied) {
    Bitboard b = 0;
    for (const auto& i : directions)
    {
        Square to = s + i.first;
        if (is_ok(to) && distance(s, to) < 4 && !(horse_path(i.first, s) & occupied))
            b |= to;
    }
    return b;
  }

}

/// safe_destination() returns the bitboard of target square for the given step
/// from the given square. If the step is off the board, returns empty bitboard.

inline Bitboard safe_destination(Square s, int step) {
    Square to = Square(s + step);
    return is_ok(to) && distance(s, to) <= 3 ? square_bb(to) : Bitboard(0);
}


/// Bitboards::pretty() returns an ASCII representation of a bitboard suitable
/// to be printed to standard output. Useful for debugging.

std::string Bitboards::pretty(Bitboard b) {

  std::string s = "+---+---+---+---+---+---+---+---+---+---+---+---+\n";

  for (Rank r = RANK_MAX; r >= RANK_1; --r)
  {
      for (File f = FILE_A; f <= FILE_MAX; ++f)
          s += b & make_square(f, r) ? "| X " : "|   ";

      s += "| " + std::to_string(1 + r) + "\n+---+---+---+---+---+---+---+---+---+---+---+---+\n";
  }
  s += "  a   b   c   d   e   f   g   h   i   j   k\n";

  return s;
}

/// Bitboards::init_pieces() initializes piece move/attack bitboards and rider types

void Bitboards::init_pieces() {

  for (PieceType pt = PAWN; pt <= KING; ++pt)
  {
      const PieceInfo* pi = pieceMap.find(pt)->second;

      // Detect rider types
      for (auto modality : {MODALITY_QUIET, MODALITY_CAPTURE})
      {
          auto& riderTypes = modality == MODALITY_CAPTURE ? AttackRiderTypes[pt] : MoveRiderTypes[pt];
          riderTypes = NO_RIDER;
          for (auto const& [d, limit] : pi->steps[modality])
          {
              if (limit && HorseDirections.find(d) != HorseDirections.end())
                  riderTypes |= RIDER_HORSE;
              if (limit && ElephantDirections.find(d) != ElephantDirections.end())
                  riderTypes |= RIDER_ELEPHANT;
          }
          for (auto const& [d, limit] : pi->slider[modality])
          {
              if (RookDirectionsH.find(d) != RookDirectionsH.end())
                  riderTypes |= RIDER_ROOK_H;
              if (RookDirectionsV.find(d) != RookDirectionsV.end())
                  riderTypes |= RIDER_ROOK_V;
          }
          for (auto const& [d, limit] : pi->hopper[modality])
          {
              if (RookDirectionsH.find(d) != RookDirectionsH.end())
                  riderTypes |= RIDER_CANNON_H;
              if (RookDirectionsV.find(d) != RookDirectionsV.end())
                  riderTypes |= RIDER_CANNON_V;
          }
      }

      // Initialize move/attack bitboards
      for (Color c : { WHITE, BLACK })
      {
          for (Square s = SQ_A1; s <= SQ_MAX; ++s)
          {
              for (auto modality : {MODALITY_QUIET, MODALITY_CAPTURE})
              {
                  auto& pseudo = modality == MODALITY_CAPTURE ? PseudoAttacks[c][pt][s] : PseudoMoves[c][pt][s];
                  auto& leaper = modality == MODALITY_CAPTURE ? LeaperAttacks[c][pt][s] : LeaperMoves[c][pt][s];
                  pseudo = 0;
                  leaper = 0;
                  for (auto const& [d, limit] : pi->steps[modality])
                  {
                      pseudo |= safe_destination(s, c == WHITE ? d : -d);
                      if (!limit)
                          leaper |= safe_destination(s, c == WHITE ? d : -d);
                  }
                  pseudo |= sliding_attack<RIDER>(pi->slider[modality], s, 0, c);
                  pseudo |= sliding_attack<UNLIMITED_RIDER>(pi->hopper[modality], s, 0, c);
              }
          }
      }


  }
}


/// Bitboards::init() initializes various bitboard tables. It is called at
/// startup and relies on global objects to be already zero-initialized.

void Bitboards::init() {

  for (unsigned i = 0; i < (1 << 16); ++i)
      PopCnt16[i] = uint8_t(std::bitset<16>(i).count());

  for (Square s = SQ_A1; s <= SQ_MAX; ++s)
      SquareBB[s] = make_bitboard(s);

  for (File f = FILE_A; f <= FILE_MAX; ++f)
      for (Rank r = RANK_1; r <= RANK_MAX; ++r)
          BoardSizeBB[f][r] = forward_file_bb(BLACK, make_square(f, r)) | SquareBB[make_square(f, r)] | (f > FILE_A ? BoardSizeBB[f - 1][r] : Bitboard(0));

  for (Square s1 = SQ_A1; s1 <= SQ_MAX; ++s1)
      for (Square s2 = SQ_A1; s2 <= SQ_MAX; ++s2)
              SquareDistance[s1][s2] = std::max(distance<File>(s1, s2), distance<Rank>(s1, s2));

#ifdef PRECOMPUTED_MAGICS
  init_magics<RIDER>(RookTableH, RookMagicsH, RookDirectionsH, RookMagicHInit);
  init_magics<RIDER>(RookTableV, RookMagicsV, RookDirectionsV, RookMagicVInit);
  init_magics<HOPPER>(CannonTableH, CannonMagicsH, RookDirectionsH, CannonMagicHInit);
  init_magics<HOPPER>(CannonTableV, CannonMagicsV, RookDirectionsV, CannonMagicVInit);
  init_magics<LAME_LEAPER>(HorseTable, HorseMagics, HorseDirections, HorseMagicInit);
  init_magics<LAME_LEAPER>(ElephantTable, ElephantMagics, ElephantDirections, ElephantMagicInit);
#else
  init_magics<RIDER>(RookTableH, RookMagicsH, RookDirectionsH);
  init_magics<RIDER>(RookTableV, RookMagicsV, RookDirectionsV);
  init_magics<HOPPER>(CannonTableH, CannonMagicsH, RookDirectionsH);
  init_magics<HOPPER>(CannonTableV, CannonMagicsV, RookDirectionsV);
  init_magics<LAME_LEAPER>(HorseTable, HorseMagics, HorseDirections);
  init_magics<LAME_LEAPER>(ElephantTable, ElephantMagics, ElephantDirections);
#endif

  init_magics<BY_HORSE>(ByHorseTable, ByHorseMagics, HorseDirections);

  init_pieces();

  for (Square s1 = SQ_A1; s1 <= SQ_MAX; ++s1)
  {
      for (PieceType pt : { ROOK })
          for (Square s2 = SQ_A1; s2 <= SQ_MAX; ++s2)
          {
              if (PseudoAttacks[WHITE][pt][s1] & s2)
              {
                  LineBB[s1][s2]    = (attacks_bb(WHITE, pt, s1, 0) & attacks_bb(WHITE, pt, s2, 0)) | s1 | s2;
                  BetweenBB[s1][s2] = (attacks_bb(WHITE, pt, s1, square_bb(s2)) & attacks_bb(WHITE, pt, s2, square_bb(s1)));
              }
              BetweenBB[s1][s2] |= s2;
          }
  }
}

namespace {

  // init_magics() computes all rook and bishop attacks at startup. Magic
  // bitboards are used to look up attacks of sliding pieces. As a reference see
  // www.chessprogramming.org/Magic_Bitboards. In particular, here we use the so
  // called "fancy" approach.

  template <MovementType MT>
  void init_magics(Bitboard table[], Magic magics[], std::map<Direction, int> directions, const Bitboard magicsInit[]) {
    Bitboard* occupancy = new Bitboard[1 << (FILE_NB + RANK_NB - 4)];
    Bitboard* reference = new Bitboard[1 << (FILE_NB + RANK_NB - 4)];
    Bitboard edges, b;
    int* epoch = new int[1 << (FILE_NB + RANK_NB - 4)]();
    int cnt = 0, size = 0;

    for (Square s = SQ_A1; s <= SQ_MAX; ++s)
    {
        // Board edges are not considered in the relevant occupancies
        edges = ((Rank1BB | rank_bb(RANK_MAX)) & ~rank_bb(s)) | ((FileABB | file_bb(FILE_MAX)) & ~file_bb(s));

        // Given a square 's', the mask is the bitboard of sliding attacks from
        // 's' computed on an empty board. The index must be big enough to contain
        // all the attacks for each possible subset of the mask and so is 2 power
        // the number of 1s of the mask. Hence we deduce the size of the shift to
        // apply to the 64 or 32 bits word to get the index.
        Magic& m = magics[s];
        // The mask for hoppers is unlimited distance, even if the hopper is limited distance (e.g., grasshopper)
        m.mask  = (MT == LAME_LEAPER ? lame_leaper_path(directions, s) : sliding_attack<MT == HOPPER ? UNLIMITED_RIDER : MT>(directions, s, 0)) & ~edges;
        
        if (HasPext)
            m.shift = popcount((m.mask << 64) >> 64);
        else
            m.shift = 128 - popcount(m.mask);

        // Set the offset for the attacks table of the square. We have individual
        // table sizes for each square with "Fancy Magic Bitboards".
        m.attacks = s == SQ_A1 ? table : magics[s - 1].attacks + size;

        // Use Carry-Rippler trick to enumerate all subsets of masks[s] and
        // store the corresponding sliding attack bitboard in reference[].
        b = size = 0;
        do {
            occupancy[size] = b;
            reference[size] = MT == LAME_LEAPER ? lame_leaper_attack(directions, s, b) : sliding_attack<MT>(directions, s, b);

            if (HasPext)
                m.attacks[pext(b, m.mask, m.shift)] = reference[size];

            size++;
            b = (b - m.mask) & m.mask;
        } while (b);

        if (HasPext)
            continue;

        // Find a magic for square 's' picking up an (almost) random number
        // until we find the one that passes the verification test.
        for (int i = 0; i < size; )
        {
            for (m.magic = 0; popcount((m.magic * m.mask) >> (SQUARE_NB - FILE_NB)) < FILE_NB - 2; )
            {
                m.magic = magicsInit[s];
            }

            // A good magic must map every possible occupancy to an index that
            // looks up the correct sliding attack in the attacks[s] database.
            // Note that we build up the database for square 's' as a side
            // effect of verifying the magic. Keep track of the attempt count
            // and save it in epoch[], little speed-up trick to avoid resetting
            // m.attacks[] after every failed attempt.
            for (++cnt, i = 0; i < size; ++i)
            {
                unsigned idx = m.index(occupancy[i]);

                if (epoch[idx] < cnt)
                {
                    epoch[idx] = cnt;
                    m.attacks[idx] = reference[i];
                }
                else if (m.attacks[idx] != reference[i])
                    break;
            }
        }
    }

    delete[] occupancy;
    delete[] reference;
    delete[] epoch;
  }

  template <MovementType MT>
  void init_magics(Bitboard table[], Magic magics[], std::map<Direction, int> directions) {

    // Optimal PRNG seeds to pick the correct magics in the shortest time
    int seeds[][RANK_NB] = { { 734, 10316, 55013, 32803, 12281, 15100,  16645, 255, 346, 89123 },
                             { 734, 10316, 55013, 32803, 12281, 15100,  16645, 255, 346, 89123 } };

    Bitboard* occupancy = new Bitboard[1 << (FILE_NB + RANK_NB - 4)];
    Bitboard* reference = new Bitboard[1 << (FILE_NB + RANK_NB - 4)];
    Bitboard edges, b;
    int* epoch = new int[1 << (FILE_NB + RANK_NB - 4)]();
    int cnt = 0, size = 0;

    for (Square s = SQ_A1; s <= SQ_MAX; ++s)
    {
        // Board edges are not considered in the relevant occupancies
        edges = ((Rank1BB | rank_bb(RANK_MAX)) & ~rank_bb(s)) | ((FileABB | file_bb(FILE_MAX)) & ~file_bb(s));

        // Given a square 's', the mask is the bitboard of sliding attacks from
        // 's' computed on an empty board. The index must be big enough to contain
        // all the attacks for each possible subset of the mask and so is 2 power
        // the number of 1s of the mask. Hence we deduce the size of the shift to
        // apply to the 64 or 32 bits word to get the index.
        Magic& m = magics[s];
        // The mask for hoppers is unlimited distance, even if the hopper is limited distance (e.g., grasshopper)
        if (MT == BY_HORSE)
            m.mask  = horse_path(directions, s);
        else
            m.mask  = (MT == LAME_LEAPER ? lame_leaper_path(directions, s) : sliding_attack<MT == HOPPER ? UNLIMITED_RIDER : MT>(directions, s, 0)) & ~edges;
        if (HasPext)
            m.shift = popcount((m.mask << 64) >> 64);
        else
            m.shift = 128 - popcount(m.mask);

        // Set the offset for the attacks table of the square. We have individual
        // table sizes for each square with "Fancy Magic Bitboards".
        m.attacks = s == SQ_A1 ? table : magics[s - 1].attacks + size;

        // Use Carry-Rippler trick to enumerate all subsets of masks[s] and
        // store the corresponding sliding attack bitboard in reference[].
        b = size = 0;
        do {
            occupancy[size] = b;
            if (MT == BY_HORSE)
                reference[size] = horse_attack(directions, s, b);
            else
                reference[size] = MT == LAME_LEAPER ? lame_leaper_attack(directions, s, b) : sliding_attack<MT>(directions, s, b);

            if (HasPext)
                m.attacks[pext(b, m.mask, m.shift)] = reference[size];

            size++;
            b = (b - m.mask) & m.mask;
        } while (b);

        if (HasPext)
            continue;

        PRNG rng(seeds[Is64Bit][rank_of(s)]);

        // Find a magic for square 's' picking up an (almost) random number
        // until we find the one that passes the verification test.
        for (int i = 0; i < size; )
        {
            for (m.magic = 0; popcount((m.magic * m.mask) >> (SQUARE_NB - FILE_NB)) < FILE_NB - 2; )
            {
                m.magic = (rng.sparse_rand<Bitboard>() << 64) ^ rng.sparse_rand<Bitboard>();
            }

            // A good magic must map every possible occupancy to an index that
            // looks up the correct sliding attack in the attacks[s] database.
            // Note that we build up the database for square 's' as a side
            // effect of verifying the magic. Keep track of the attempt count
            // and save it in epoch[], little speed-up trick to avoid resetting
            // m.attacks[] after every failed attempt.
            for (++cnt, i = 0; i < size; ++i)
            {
                unsigned idx = m.index(occupancy[i]);

                if (epoch[idx] < cnt)
                {
                    epoch[idx] = cnt;
                    m.attacks[idx] = reference[i];
                }
                else if (m.attacks[idx] != reference[i])
                    break;
            }
        }
    }

    delete[] occupancy;
    delete[] reference;
    delete[] epoch;
  }
}

} // namespace Stockfish

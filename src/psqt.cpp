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


#include "psqt.h"

#include <algorithm>
#include <sstream>
#include <math.h>

#include "bitboard.h"
#include "types.h"

#include "piece.h"
#include "variant.h"
#include "misc.h"


namespace Stockfish {

Value EvalPieceValue[PHASE_NB][PIECE_NB];
Value CapturePieceValue[PHASE_NB][PIECE_NB];

Value PieceValue[PHASE_NB][PIECE_NB] = {
  {
    VALUE_ZERO, PawnValueMg, KnightValueMg, BishopValueMg, RookValueMg, QueenValueMg, FersValueMg, AlfilValueMg,
    FersAlfilValueMg, SilverValueMg, AiwokValueMg, BersValueMg, ArchbishopValueMg, ChancellorValueMg, AmazonValueMg, KnibisValueMg,
    BiskniValueMg, KnirooValueMg, RookniValueMg, ShogiPawnValueMg, LanceValueMg, ShogiKnightValueMg, GoldValueMg, DragonHorseValueMg,
    ClobberPieceValueMg, BreakthroughPieceValueMg, ImmobilePieceValueMg, CannonPieceValueMg, JanggiCannonPieceValueMg, SoldierValueMg, HorseValueMg, ElephantValueMg,
    JanggiElephantValueMg, BannerValueMg, WazirValueMg, CommonerValueMg, CentaurValueMg, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO,
    VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO,
    VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO,
    VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO,

    VALUE_ZERO, PawnValueMg, KnightValueMg, BishopValueMg, RookValueMg, QueenValueMg, FersValueMg, AlfilValueMg,
    FersAlfilValueMg, SilverValueMg, AiwokValueMg, BersValueMg, ArchbishopValueMg, ChancellorValueMg, AmazonValueMg, KnibisValueMg,
    BiskniValueMg, KnirooValueMg, RookniValueMg, ShogiPawnValueMg, LanceValueMg, ShogiKnightValueMg, GoldValueMg, DragonHorseValueMg,
    ClobberPieceValueMg, BreakthroughPieceValueMg, ImmobilePieceValueMg, CannonPieceValueMg, JanggiCannonPieceValueMg, SoldierValueMg, HorseValueMg, ElephantValueMg,
    JanggiElephantValueMg, BannerValueMg, WazirValueMg, CommonerValueMg, CentaurValueMg, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO,
    VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO,
    VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO,
    VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO,
  },
  {
    VALUE_ZERO, PawnValueEg, KnightValueEg, BishopValueEg, RookValueEg, QueenValueEg, FersValueEg, AlfilValueEg,
    FersAlfilValueEg, SilverValueEg, AiwokValueEg, BersValueEg, ArchbishopValueEg, ChancellorValueEg, AmazonValueEg, KnibisValueEg,
    BiskniValueEg, KnirooValueEg, RookniValueEg, ShogiPawnValueEg, LanceValueEg, ShogiKnightValueEg, GoldValueEg, DragonHorseValueEg,
    ClobberPieceValueEg, BreakthroughPieceValueEg, ImmobilePieceValueEg, CannonPieceValueEg, JanggiCannonPieceValueEg, SoldierValueEg, HorseValueEg, ElephantValueEg,
    JanggiElephantValueEg, BannerValueEg, WazirValueEg, CommonerValueEg, CentaurValueEg, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO,
    VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO,
    VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO,
    VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO,

    VALUE_ZERO, PawnValueEg, KnightValueEg, BishopValueEg, RookValueEg, QueenValueEg, FersValueEg, AlfilValueEg,
    FersAlfilValueEg, SilverValueEg, AiwokValueEg, BersValueEg, ArchbishopValueEg, ChancellorValueEg, AmazonValueEg, KnibisValueEg,
    BiskniValueEg, KnirooValueEg, RookniValueEg, ShogiPawnValueEg, LanceValueEg, ShogiKnightValueEg, GoldValueEg, DragonHorseValueEg,
    ClobberPieceValueEg, BreakthroughPieceValueEg, ImmobilePieceValueEg, CannonPieceValueEg, JanggiCannonPieceValueEg, SoldierValueEg, HorseValueEg, ElephantValueEg,
    JanggiElephantValueEg, BannerValueEg, WazirValueEg, CommonerValueEg, CentaurValueEg, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO,
    VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO,
    VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO,
    VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO, VALUE_ZERO,
  },
};


namespace
{

auto constexpr S = make_score;

// 'Bonus' contains Piece-Square parameters.
// Scores are explicit for files A to D, implicitly mirrored for E to H.
constexpr Score Bonus[PIECE_TYPE_NB][RANK_NB][int(FILE_NB) / 2] = {
  { },
  { },
  { // Knight
   { S(-175, -96), S(-92,-65), S(-74,-49), S(-73,-21) },
   { S( -77, -67), S(-41,-54), S(-27,-18), S(-15,  8) },
   { S( -61, -40), S(-17,-27), S(  6, -8), S( 12, 29) },
   { S( -35, -35), S(  8, -2), S( 40, 13), S( 49, 28) },
   { S( -34, -45), S( 13,-16), S( 44,  9), S( 51, 39) },
   { S(  -9, -51), S( 22,-44), S( 58,-16), S( 53, 17) },
   { S( -67, -69), S(-27,-50), S(  4,-51), S( 37, 12) },
   { S(-201,-100), S(-83,-88), S(-56,-56), S(-26,-17) }
  },
  { // Bishop
   { S(-37,-40), S(-4 ,-21), S( -6,-26), S(-16, -8) },
   { S(-11,-26), S(  6, -9), S( 13,-12), S(  3,  1) },
   { S(-5 ,-11), S( 15, -1), S( -4, -1), S( 12,  7) },
   { S(-4 ,-14), S(  8, -4), S( 18,  0), S( 27, 12) },
   { S(-8 ,-12), S( 20, -1), S( 15,-10), S( 22, 11) },
   { S(-11,-21), S(  4,  4), S(  1,  3), S(  8,  4) },
   { S(-12,-22), S(-10,-14), S(  4, -1), S(  0,  1) },
   { S(-34,-32), S(  1,-29), S(-10,-26), S(-16,-17) }
  },
  { // Rook
   { S(-31, -9), S(-20,-13), S(-14,-10), S(-5, -9) },
   { S(-21,-12), S(-13, -9), S( -8, -1), S( 6, -2) },
   { S(-25,  6), S(-11, -8), S( -1, -2), S( 3, -6) },
   { S(-13, -6), S( -5,  1), S( -4, -9), S(-6,  7) },
   { S(-27, -5), S(-15,  8), S( -4,  7), S( 3, -6) },
   { S(-22,  6), S( -2,  1), S(  6, -7), S(12, 10) },
   { S( -2,  4), S( 12,  5), S( 16, 20), S(18, -5) },
   { S(-17, 18), S(-19,  0), S( -1, 19), S( 9, 13) }
  },
  { // Queen
   { S( 3,-69), S(-5,-57), S(-5,-47), S( 4,-26) },
   { S(-3,-54), S( 5,-31), S( 8,-22), S(12, -4) },
   { S(-3,-39), S( 6,-18), S(13, -9), S( 7,  3) },
   { S( 4,-23), S( 5, -3), S( 9, 13), S( 8, 24) },
   { S( 0,-29), S(14, -6), S(12,  9), S( 5, 21) },
   { S(-4,-38), S(10,-18), S( 6,-11), S( 8,  1) },
   { S(-5,-50), S( 6,-27), S(10,-24), S( 8, -8) },
   { S(-2,-74), S(-2,-52), S( 1,-43), S(-2,-34) }
  }
};

constexpr Score KingBonus[RANK_NB][int(FILE_NB) / 2] = {
   { S(271,  1), S(327, 45), S(271, 85), S(198, 76) },
   { S(278, 53), S(303,100), S(234,133), S(179,135) },
   { S(195, 88), S(258,130), S(169,169), S(120,175) },
   { S(164,103), S(190,156), S(138,172), S( 98,172) },
   { S(154, 96), S(179,166), S(105,199), S( 70,199) },
   { S(123, 92), S(145,172), S( 81,184), S( 31,191) },
   { S( 88, 47), S(120,121), S( 65,116), S( 33,131) },
   { S( 59, 11), S( 89, 59), S( 45, 73), S( -1, 78) }
};

constexpr Score PBonus[RANK_NB][FILE_NB] =
  { // Pawn (asymmetric distribution)
   { },
   { S(  2, -8), S(  4, -6), S( 11,  9), S( 18,  5), S( 16, 16), S( 21,  6), S(  9, -6), S( -3,-18) },
   { S( -9, -9), S(-15, -7), S( 11,-10), S( 15,  5), S( 31,  2), S( 23,  3), S(  6, -8), S(-20, -5) },
   { S( -3,  7), S(-20,  1), S(  8, -8), S( 19, -2), S( 39,-14), S( 17,-13), S(  2,-11), S( -5, -6) },
   { S( 11, 12), S( -4,  6), S(-11,  2), S(  2, -6), S( 11, -5), S(  0, -4), S(-12, 14), S(  5,  9) },
   { S(  3, 27), S(-11, 18), S( -6, 19), S( 22, 29), S( -8, 30), S( -5,  9), S(-14,  8), S(-11, 14) },
   { S( -7, -1), S(  6,-14), S( -2, 13), S(-11, 22), S(  4, 24), S(-14, 17), S( 10,  7), S( -9,  7) }
  };

// Estimate piece value
Value piece_value(Phase phase, PieceType pt)
{
    const PieceInfo* pi = pieceMap.find(pt)->second;
    int v0 =  (phase == MG ?  55 :  60) * pi->steps[MODALITY_CAPTURE].size()
            + (phase == MG ?  30 :  40) * pi->steps[MODALITY_QUIET].size()
            + (phase == MG ? 185 : 180) * pi->slider[MODALITY_CAPTURE].size()
            + (phase == MG ?  55 :  50) * pi->slider[MODALITY_QUIET].size()
            // Hoppers are more useful with more pieces on the board
            + (phase == MG ? 100 :  80) * pi->hopper[MODALITY_CAPTURE].size()
            + (phase == MG ?  80 :  60) * pi->hopper[MODALITY_QUIET].size()
            // Rook sliding directions are more valuable, especially in endgame
            + (phase == MG ?  10 :  30) * std::count_if(pi->slider[MODALITY_CAPTURE].begin(), pi->slider[MODALITY_CAPTURE].end(), [](const std::pair<const Direction, int>& d) { return std::abs(d.first) == NORTH || std::abs(d.first) == 1; })
            + (phase == MG ?  30 :  45) * std::count_if(pi->slider[MODALITY_QUIET].begin(), pi->slider[MODALITY_QUIET].end(), [](const std::pair<const Direction, int>& d) { return std::abs(d.first) == NORTH || std::abs(d.first) == 1; });
    return Value(v0 * exp(double(v0) / 10000));
}

} // namespace


namespace PSQT
{

Score psq[PIECE_NB][SQUARE_NB + 1];

// PSQT::init() initializes piece-square tables: the white halves of the tables are
// copied from Bonus[] and PBonus[], adding the piece value, then the black halves of
// the tables are initialized by flipping and changing the sign of the white scores.
void init(const Variant* v) {

  PieceType strongestPiece = NO_PIECE_TYPE;
  for (PieceType pt : pieceTypes)
  {
      if (is_custom(pt))
      {
          PieceValue[MG][pt] = piece_value(MG, pt);
          PieceValue[EG][pt] = piece_value(EG, pt);
      }

      if (PieceValue[MG][pt] > PieceValue[MG][strongestPiece])
          strongestPiece = pt;
  }

  for (PieceType pt = PAWN; pt <= KING; ++pt)
  {
      Piece pc = make_piece(WHITE, pt);

      Score score = make_score(PieceValue[MG][pc], PieceValue[EG][pc]);
      
      const PieceInfo* pi = pieceMap.find(pt)->second;
      bool isSlider = pi->slider[MODALITY_QUIET].size() || pi->slider[MODALITY_CAPTURE].size() || pi->hopper[MODALITY_QUIET].size() || pi->hopper[MODALITY_CAPTURE].size();
      bool isPawn = !isSlider && pi->steps[MODALITY_QUIET].size() && !std::any_of(pi->steps[MODALITY_QUIET].begin(), pi->steps[MODALITY_QUIET].end(), [](const std::pair<const Direction, int>& d) { return d.first < SOUTH / 2; });
      bool isSlowLeaper = !isSlider && !std::any_of(pi->steps[MODALITY_QUIET].begin(), pi->steps[MODALITY_QUIET].end(), [](const std::pair<const Direction, int>& d) { return dist(d.first) > 1; });

      // Scale slider piece values with board size
      if (isSlider)
      {
          constexpr int lc = 5;
          constexpr int rm = 5;
          constexpr int r0 = rm + RANK_8;
          int r1 = rm + (v->maxRank + v->maxFile) / 2;
          int leaper = pi->steps[MODALITY_QUIET].size() + pi->steps[MODALITY_CAPTURE].size();
          int slider = pi->slider[MODALITY_QUIET].size() + pi->slider[MODALITY_CAPTURE].size() + pi->hopper[MODALITY_QUIET].size() + pi->hopper[MODALITY_CAPTURE].size();
          score = make_score(mg_value(score) * (lc * leaper + r1 * slider) / (lc * leaper + r0 * slider),
                             eg_value(score) * (lc * leaper + r1 * slider) / (lc * leaper + r0 * slider));
      }

      CapturePieceValue[MG][pc] = CapturePieceValue[MG][~pc] = mg_value(score);
      CapturePieceValue[EG][pc] = CapturePieceValue[EG][~pc] = eg_value(score);

      EvalPieceValue[MG][pc] = EvalPieceValue[MG][~pc] = mg_value(score);
      EvalPieceValue[EG][pc] = EvalPieceValue[EG][~pc] = eg_value(score);

      // Determine pawn rank
      std::istringstream ss(v->startFen);
      unsigned char token;
      Rank rc = v->maxRank;
      Rank pawnRank = RANK_2;
      while ((ss >> token) && !isspace(token))
      {
          if (token == '/')
              --rc;
          else if (token == v->pieceToChar[PAWN])
              pawnRank = rc;
      }

      for (Square s = SQ_A1; s <= SQ_MAX; ++s)
      {
          File f = std::max(File(edge_distance(file_of(s), v->maxFile)), FILE_A);
          Rank r = rank_of(s);
          psq[ pc][s] = score + (  pt == KING  ? KingBonus[std::clamp(Rank(r - pawnRank + 1), RANK_1, RANK_8)][std::min(f, FILE_D)] * 1
                                 : pt <= QUEEN ? Bonus[pc][std::min(r, RANK_8)][std::min(f, FILE_D)]
                                 : pt == HORSE ? Bonus[KNIGHT][std::min(r, RANK_8)][std::min(f, FILE_D)]
                                 : isSlider    ? make_score(5, 5) * (2 * f + std::max(std::min(r, Rank(v->maxRank - r)), RANK_1) - v->maxFile - 1)
                                 : isPawn      ? make_score(5, 5) * (2 * f - v->maxFile)
                                               : make_score(10, 10) * (1 + isSlowLeaper) * (f + std::max(std::min(r, Rank(v->maxRank - r)), RANK_1) - v->maxFile / 2));
          // Add a penalty for unpromoted soldiers
          if (pt == SOLDIER && r < RANK_6)
              psq[pc][s] -= score * (RANK_6 - r) / (4 + f);
          psq[~pc][rank_of(s) <= v->maxRank ? flip_rank(s, v->maxRank) : s] = -psq[pc][s];
      }
      // Pieces in hand
      psq[ pc][SQ_NONE] = score + make_score(35, 10) * (1 + !isSlider);
      psq[~pc][SQ_NONE] = -psq[pc][SQ_NONE];
  }
}

} // namespace PSQT

} // namespace Stockfish

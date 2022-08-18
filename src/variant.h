/*
  Fairy-Stockfish, a UCI chess variant playing engine derived from Stockfish
  Copyright (C) 2018-2022 Fabian Fichter

  Fairy-Stockfish is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  Fairy-Stockfish is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef VARIANT_H_INCLUDED
#define VARIANT_H_INCLUDED

#include <array>
#include <set>
#include <map>
#include <vector>
#include <string>
#include <functional>
#include <sstream>
#include <iostream>

#include "types.h"
#include "bitboard.h"

namespace Stockfish {

/// Variant struct stores information needed to determine the rules of a variant.

constexpr std::array<PieceType, 7> pieceTypes = { ROOK, FERS, CANNON, SOLDIER, HORSE, ELEPHANT, KING };

struct Variant {
  std::string variantTemplate = "fairy";
  std::string pieceToCharTable = "-";
  Rank maxRank = RANK_10;
  File maxFile = FILE_I;
  std::string pieceToChar = std::string(PIECE_NB, ' ');
  std::string pieceToCharSynonyms = std::string(PIECE_NB, ' ');
  std::string startFen = "rnbakabnr/9/1c5c1/p1p1p1p1p/9/9/P1P1P1P1P/1C5C1/9/RNBAKABNR w - - 0 1";
  Bitboard mobilityRegion[COLOR_NB][PIECE_TYPE_NB] = {};

  std::string nnueAlias = "";
  int nnueDimensions;
  int pieceSquareIndex[COLOR_NB][PIECE_NB];
  int pieceHandIndex[COLOR_NB][PIECE_NB];
  int kingSquareIndex[SQUARE_NB];
  int nnueMaxPieces;

  void add_piece(PieceType pt, char c, char c2 = ' ') {
      pieceToChar[make_piece(WHITE, pt)] = toupper(c);
      pieceToChar[make_piece(BLACK, pt)] = tolower(c);
      pieceToCharSynonyms[make_piece(WHITE, pt)] = toupper(c2);
      pieceToCharSynonyms[make_piece(BLACK, pt)] = tolower(c2);
  }

  // Reset values that always need to be redefined
  Variant* init() {
      nnueAlias = "";
      return this;
  }

  // Pre-calculate derived properties
  Variant* conclude() {

      // Initialize calculated NNUE properties
      int nnueSquares = 90;
      int nnuePockets = 0;
      int nnueNonDropPieceIndices = (2 * pieceTypes.size() - (KING != NO_PIECE_TYPE)) * nnueSquares;
      int nnuePieceIndices = nnueNonDropPieceIndices + 2 * (pieceTypes.size() - (KING != NO_PIECE_TYPE)) * nnuePockets;
      int i = 0;
      for (PieceType pt : pieceTypes)
      {
          for (Color c : { WHITE, BLACK})
          {
              pieceSquareIndex[c][make_piece(c, pt)] = 2 * i * nnueSquares;
              pieceSquareIndex[c][make_piece(~c, pt)] = (2 * i + (pt != KING)) * nnueSquares;
              pieceHandIndex[c][make_piece(c, pt)] = 2 * i * nnuePockets + nnueNonDropPieceIndices;
              pieceHandIndex[c][make_piece(~c, pt)] = (2 * i + 1) * nnuePockets + nnueNonDropPieceIndices;
          }
          i++;
      }

      // Map king squares to enumeration of actually available squares.
      // E.g., for xiangqi map from 0-89 to 0-8.
      // Variants might be initialized before bitboards, so do not rely on precomputed bitboards (like SquareBB).
      int nnueKingSquare = 0;
      if (KING)
          for (Square s = SQ_A1; s < nnueSquares; ++s)
          {
              Square bitboardSquare = Square(s + s / (maxFile + 1) * (FILE_MAX - maxFile));
              if (   !mobilityRegion[WHITE][KING] || !mobilityRegion[BLACK][KING]
                  || (mobilityRegion[WHITE][KING] & make_bitboard(bitboardSquare))
                  || (mobilityRegion[BLACK][KING] & make_bitboard(relative_square(BLACK, bitboardSquare, RANK_10))))
              {
                  kingSquareIndex[s] = nnueKingSquare++ * nnuePieceIndices;
              }
          }
      else
          kingSquareIndex[SQ_A1] = nnueKingSquare++ * nnuePieceIndices;
      nnueDimensions = nnueKingSquare * nnuePieceIndices;

      // Determine maximum piece count
      std::istringstream ss(startFen);
      ss >> std::noskipws;
      unsigned char token;
      nnueMaxPieces = 0;
      while ((ss >> token) && !isspace(token))
      {
          if (pieceToChar.find(token) != std::string::npos || pieceToCharSynonyms.find(token) != std::string::npos)
              nnueMaxPieces++;
      }

      return this;
  }

  void late_init() {
      for (Color c : {WHITE, BLACK})
          for (int pt = 0; pt < PIECE_TYPE_NB; pt++) {
              auto board_bb = board_size_bb(maxFile, maxRank);
              BoardBB[c][pt] = mobilityRegion[c][pt] ? mobilityRegion[c][pt] & board_bb : board_bb;
          }
  }
};

class VariantMap : public std::map<std::string, const Variant*> {
public:
  void init();
  void clear_all();
  std::vector<std::string> get_keys();

private:
  void add(std::string s, Variant* v);
};

extern VariantMap variants;

} // namespace Stockfish

#endif // #ifndef VARIANT_H_INCLUDED

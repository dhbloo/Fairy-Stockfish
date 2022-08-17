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

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include "piece.h"
#include "variant.h"

using std::string;

namespace Stockfish {

VariantMap variants; // Global object

namespace {
    Variant* xiangqi_variant() {
        Variant* v = new Variant();
        v->variantTemplate = "xiangqi";
        v->pieceToCharTable = "PN.R.AB..K.C..........pn.r.ab..k.c..........";
        v->reset_pieces();
        v->add_piece(ROOK, 'r');
        v->add_piece(HORSE, 'n', 'h');
        v->add_piece(KING, 'k');
        v->add_piece(CANNON, 'c');
        v->add_piece(SOLDIER, 'p');
        v->add_piece(ELEPHANT, 'b', 'e');
        v->add_piece(FERS, 'a');
        v->mobilityRegion[WHITE][KING] = (Rank1BB | Rank2BB | Rank3BB) & (FileDBB | FileEBB | FileFBB);
        v->mobilityRegion[BLACK][KING] = (Rank8BB | Rank9BB | Rank10BB) & (FileDBB | FileEBB | FileFBB);
        v->mobilityRegion[WHITE][FERS] = v->mobilityRegion[WHITE][KING];
        v->mobilityRegion[BLACK][FERS] = v->mobilityRegion[BLACK][KING];
        v->mobilityRegion[WHITE][ELEPHANT] = Rank1BB | Rank2BB | Rank3BB | Rank4BB | Rank5BB;
        v->mobilityRegion[BLACK][ELEPHANT] = Rank6BB | Rank7BB | Rank8BB | Rank9BB | Rank10BB;
        v->kingType = WAZIR;
        v->promotionPieceTypes = {};
        v->stalemateValue = -VALUE_MATE;
        v->soldierPromotionRank = RANK_6;
        return v;
    }

} // namespace


/// VariantMap::init() is called at startup to initialize all predefined variants

void VariantMap::init() {
    add("xiangqi", xiangqi_variant());
}

void VariantMap::add(std::string s, Variant* v) {
  insert(std::pair<std::string, const Variant*>(s, v->conclude()));
}

void VariantMap::clear_all() {
  for (auto const& element : *this)
      delete element.second;
  clear();
}

std::vector<std::string> VariantMap::get_keys() {
  std::vector<std::string> keys;
  for (auto const& element : *this)
      keys.push_back(element.first);
  return keys;
}

} // namespace Stockfish

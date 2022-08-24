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
#include <cstdlib>
#include <cstring>   // For std::memset
#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <streambuf>
#include <vector>

#include "bitboard.h"
#include "evaluate.h"
#include "misc.h"
#include "thread.h"
#include "timeman.h"
#include "uci.h"
#include "incbin/incbin.h"


// Macro to embed the default efficiently updatable neural network (NNUE) file
// data in the engine binary (using incbin.h, by Dale Weiler).
// This macro invocation will declare the following three variables
//     const unsigned char        gEmbeddedNNUEData[];  // a pointer to the embedded data
//     const unsigned char *const gEmbeddedNNUEEnd;     // a marker to the end
//     const unsigned int         gEmbeddedNNUESize;    // the size of the embedded file
// Note that this does not work in Microsoft Visual Studio.
#if !defined(_MSC_VER) && !defined(NNUE_EMBEDDING_OFF)
  INCBIN(EmbeddedNNUE, EvalFileDefaultName);
#else
  const unsigned char        gEmbeddedNNUEData[1] = {0x0};
  [[maybe_unused]]
  const unsigned char *const gEmbeddedNNUEEnd = &gEmbeddedNNUEData[1];
  const unsigned int         gEmbeddedNNUESize = 1;
#endif


using namespace std;

namespace Stockfish {

const Variant* currentNnueVariant;

namespace Eval {

  string eval_file_loaded = "None";

  /// NNUE::init() tries to load a NNUE network at startup time, or when the engine
  /// receives a UCI command "setoption name EvalFile value nn-[a-z0-9]{12}.nnue"
  /// The name of the NNUE network is always retrieved from the EvalFile option.
  /// We search the given network in three locations: internally (the default
  /// network may be embedded in the binary), in the active working directory and
  /// in the engine directory. Distro packagers may define the DEFAULT_NNUE_DIRECTORY
  /// variable to have the engine search in a special directory in their distro.

  void NNUE::init() {

    string eval_file = string(Options["EvalFile"]);

    // Restrict NNUE usage to corresponding variant
    // Support multiple variant networks separated by semicolon(Windows)/colon(Unix)
    stringstream ss(eval_file);
    string variant = string("xiangqi");
    bool foundNNUE = false;
    while (getline(ss, eval_file, UCI::SepChar))
    {
        string basename = eval_file.substr(eval_file.find_last_of("\\/") + 1);
        string nnueAlias = variants.find(variant)->second->nnueAlias;
        if (basename.rfind(variant, 0) != string::npos || (!nnueAlias.empty() && basename.rfind(nnueAlias, 0) != string::npos))
        {
            foundNNUE = true;
            break;
        }
    }
    if (!foundNNUE)
        return;

    currentNnueVariant = variants.find(variant)->second;

    #if defined(DEFAULT_NNUE_DIRECTORY)
    #define stringify2(x) #x
    #define stringify(x) stringify2(x)
    vector<string> dirs = { "<internal>" , "" , CommandLine::binaryDirectory , stringify(DEFAULT_NNUE_DIRECTORY) };
    #else
    vector<string> dirs = { "<internal>" , "" , CommandLine::binaryDirectory };
    #endif

    for (string directory : dirs)
        if (eval_file_loaded != eval_file)
        {
            if (directory != "<internal>")
            {
                ifstream stream(directory + eval_file, ios::binary);
                if (load_eval(eval_file, stream))
                    eval_file_loaded = eval_file;
            }

            if (directory == "<internal>" && eval_file == EvalFileDefaultName)
            {
                // C++ way to prepare a buffer for a memory stream
                class MemoryBuffer : public basic_streambuf<char> {
                    public: MemoryBuffer(char* p, size_t n) { setg(p, p, p + n); setp(p, p + n); }
                };

                MemoryBuffer buffer(const_cast<char*>(reinterpret_cast<const char*>(gEmbeddedNNUEData)),
                                    size_t(gEmbeddedNNUESize));

                istream stream(&buffer);
                if (load_eval(eval_file, stream))
                    eval_file_loaded = eval_file;
            }
        }
  }

  /// NNUE::verify() verifies that the last net used was loaded successfully
  void NNUE::verify() {

    string eval_file = string(Options["EvalFile"]);

    if (eval_file.find(eval_file_loaded) == string::npos)
    {
        UCI::OptionsMap defaults;
        UCI::init(defaults);

        string msg1 = "Network evaluation parameters compatible with the engine must be available.";
        string msg2 = "The network file " + eval_file + " was not loaded successfully.";
        string msg3 = "The UCI option EvalFile might need to specify the full path, including the directory name, to the network file.";
        string msg4 = "The default net can be downloaded from QQ group: 755655813.";
        string msg5 = "The engine will be terminated now.";

        sync_cout << "info string ERROR: " << msg1 << sync_endl;
        sync_cout << "info string ERROR: " << msg2 << sync_endl;
        sync_cout << "info string ERROR: " << msg3 << sync_endl;
        sync_cout << "info string ERROR: " << msg4 << sync_endl;
        sync_cout << "info string ERROR: " << msg5 << sync_endl;

        exit(EXIT_FAILURE);
    }

    sync_cout << "info string NNUE evaluation using " << eval_file_loaded << " enabled" << sync_endl;
  }
}

namespace Trace {

  enum Tracing { NO_TRACE, TRACE };

  enum Term { // The first PIECE_TYPE_NB entries are reserved for PieceType
    MATERIAL = PIECE_TYPE_NB, IMBALANCE, MOBILITY, THREAT, PASSED, SPACE, VARIANT, WINNABLE, TOTAL, TERM_NB
  };

  Score scores[TERM_NB][COLOR_NB];

  double to_cp(Value v) { return double(v) / PawnValueEg; }

  void add(int idx, Color c, Score s) {
    scores[idx][c] = s;
  }

  void add(int idx, Score w, Score b = SCORE_ZERO) {
    scores[idx][WHITE] = w;
    scores[idx][BLACK] = b;
  }

  std::ostream& operator<<(std::ostream& os, Score s) {
    os << std::setw(5) << to_cp(mg_value(s)) << " "
       << std::setw(5) << to_cp(eg_value(s));
    return os;
  }

  std::ostream& operator<<(std::ostream& os, Term t) {

    if (t == MATERIAL || t == IMBALANCE || t == WINNABLE || t == TOTAL)
        os << " ----  ----"    << " | " << " ----  ----";
    else
        os << scores[t][WHITE] << " | " << scores[t][BLACK];

    os << " | " << scores[t][WHITE] - scores[t][BLACK] << " |\n";
    return os;
  }
}

using namespace Trace;

/// evaluate() is the evaluator for the outer world. It returns a static
/// evaluation of the position from the point of view of the side to move.

Value Eval::evaluate(const Position& pos, int* complexity) {

  Color stm = pos.side_to_move();
  Value psq = pos.psq_eg_stm();

  int nnueComplexity;
  int scale = 1037 + 53 * pos.non_pawn_material() / 4096;
  Value optimism = pos.this_thread()->optimism[stm];

  Value nnue = NNUE::evaluate(pos, &nnueComplexity);
  // Blend nnue complexity with (semi)classical complexity
  nnueComplexity = (110 * nnueComplexity + 148 * abs(nnue - psq)) / 256;
  if (complexity) // Return hybrid NNUE complexity to caller
      *complexity = nnueComplexity;

  optimism = optimism * (239 + nnueComplexity) / 256;
  Value v = (nnue * scale + optimism * (scale - 754)) / 1024;

  // Guarantee evaluation does not hit the tablebase range
  v = std::clamp(v, VALUE_TB_LOSS_IN_MAX_PLY + 1, VALUE_TB_WIN_IN_MAX_PLY - 1);

  return v;
}


/// trace() is like evaluate(), but instead of returning a value, it returns
/// a string (suitable for outputting to stdout) that contains the detailed
/// descriptions and values of each evaluation term. Useful for debugging.
/// Trace scores are from white's point of view

std::string Eval::trace(Position& pos) {

  if (pos.checkers())
      return "Final evaluation: none (in check)";

  std::stringstream ss;
  ss << std::showpoint << std::noshowpos << std::fixed << std::setprecision(2);

  Value v;

  std::memset(scores, 0, sizeof(scores));

  pos.this_thread()->trend = SCORE_ZERO; // Reset any dynamic contempt

  ss << '\n' << NNUE::trace(pos) << '\n';

  ss << std::showpoint << std::showpos << std::fixed << std::setprecision(2) << std::setw(15);
  v = NNUE::evaluate(pos);
  v = pos.side_to_move() == WHITE ? v : -v;
  ss << "NNUE evaluation        " << to_cp(v) << " (white side)\n";

  v = evaluate(pos);
  v = pos.side_to_move() == WHITE ? v : -v;
  ss << "Final evaluation       " << to_cp(v) << " (white side)";
  if (pos.nnue_applicable())
     ss << " [with scaled NNUE, hybrid, ...]";
  ss << "\n";

  return ss.str();
}

} // namespace Stockfish

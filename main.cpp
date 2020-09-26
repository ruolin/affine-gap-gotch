/*
The 1982 Gotoh algorithm has flaws. The implementation address this flaws. 
Inspired by paper https://www.biorxiv.org/content/10.1101/031500v1.full.pdf
and the correct algorithm is documented at 
https://www.researchgate.net/publication/19580571_Optimal_sequence_alignment_using_affine_gap_costs
*/

#include <iostream>
#include <string>
#include <vector>
#include <climits>
#include <array>
#include <assert.h>

//#include "util.h"

using std::string;
using std::vector;

class AlignPath {
  using Node = std::pair<int,int>;
  //std::vector<Node> nodes_;
  const std::string ref_;
  const std::string alt_;
  Node head_;
  Node previous_;
  std::array<std::string, 2> align_;
public:
  AlignPath(const std::string& r, const std::string& a): ref_(r), alt_(a),
            previous_ (std::pair<int,int>(alt_.length() + 1, ref_.length() + 1)),
            head_(std::pair<int,int>(alt_.length(), ref_.length())) {}

  void GoNext(int i, int j) {
    assert (i >=0 && j >=0);
    //std::cerr << "i, j " << i<<", " <<j <<std::endl;
    //std::cerr << "head " << head_.first<<", " <<head_.second <<std::endl;
    if (i == head_.first && j + 1== head_.second ) {
      align_[0] += ref_[j];
      align_[1] += '-';

    } else if (i + 1 == head_.first && j == head_.second) {
      align_[0] += '-';
      align_[1] += alt_[i];

    } else if (i + 1 == head_.first && j + 1== head_.second) {
      align_[0] += ref_[j];
      align_[1] += alt_[i];
    }
    else {
      assert (false);
    }
    previous_ = head_;
    head_ = std::pair<int,int>(i,j);
  }

  const std::string& Ref() const {return align_[0];}
  const std::string& Alt() const {return align_[1];}
};

class AffineGap {
  /*
   * matrix as Reference at row top Query at column left
   */
  string query_;
  string ref_;
  AlignPath align_path_;
  int nrow_;
  int ncol_;
  vector<vector<int>> R_; // match matrix
  vector<vector<int>> P_; // vertical insertion matrix
  vector<vector<int>> Q_; // horizontal deletion matrix
  vector<vector<bool>> vert_whole_; // a
  vector<vector<bool>> hori_whole_; // b
  vector<vector<bool>> diag_whole_; // c
  vector<vector<bool>> vert_top_half_; // d
  vector<vector<bool>> vert_bottom_half_; // e
  vector<vector<bool>> hori_left_half_; // f
  vector<vector<bool>> hori_right_half_; // g

  //gap score = gap_open + gap_ext * gap_len
  const static int gap_open_ = -5;
  const static int gap_ext_ = -1;
  static int DiagScore_(const char& a, const char& b) {
    return a == b? 0: -1;
  }

public:
  AffineGap(const string& ref, const string& query): query_(query), ref_(ref), align_path_(ref_, query_),
                                                     nrow_((int) query_.length() + 1),
                                                     ncol_((int) ref_.length() + 1),
                                                     R_(vector<vector<int>>(nrow_, vector<int>(ncol_))),
                                                     P_(vector<vector<int>>(nrow_, vector<int>(ncol_))),
                                                     Q_(vector<vector<int>>(nrow_, vector<int>(ncol_))),
                                                     vert_whole_(vector<vector<bool>>(nrow_ + 1, vector<bool>(ncol_ + 1))),
                                                     hori_whole_(vector<vector<bool>>(nrow_ + 1, vector<bool>(ncol_ + 1))),
                                                     diag_whole_(vector<vector<bool>>(nrow_ + 1, vector<bool>(ncol_ + 1))),
                                                     vert_top_half_(vector<vector<bool>>(nrow_ + 1, vector<bool>(ncol_ + 1))),
                                                     vert_bottom_half_(vector<vector<bool>>(nrow_ + 1, vector<bool>(ncol_ + 1))),
                                                     hori_left_half_(vector<vector<bool>>(nrow_ + 1, vector<bool>(ncol_ + 1))),
                                                     hori_right_half_(vector<vector<bool>>(nrow_ + 1, vector<bool>(ncol_ + 1)))
  {
    // init
    for (int j = 0; j < ncol_; ++j) {
      P_[0][j] = 2 * gap_open_ + std::max(ncol_, nrow_) * gap_ext_ - 1; // ensure a large number
      R_[0][j] = gap_open_ + j * gap_ext_;
    }
    for (int i = 0; i < nrow_; ++i) {
      Q_[i][0] = 2 * gap_open_ + std::max(ncol_, nrow_) * gap_ext_ - 1; // ensure a large number
      R_[i][0] = gap_open_ + i * gap_ext_;
    }
    R_[0][0] = 0;
    diag_whole_[nrow_][ncol_] = 1;
   //
    for (int i = 0; i < nrow_; ++i) {
      for (int j = 0; j < ncol_; ++j) {
        if (i != 0) {
          P_[i][j] = gap_ext_ + std::max(P_[i-1][j], R_[i-1][j] + gap_open_);
          if (P_[i][j] == gap_ext_ + P_[i-1][j]) vert_top_half_[i-1][j] = 1;
          if (P_[i][j] == gap_ext_ + gap_open_ + R_[i-1][j]) vert_bottom_half_[i-1][j] = 1;
        }
        if (j != 0) {
          Q_[i][j] = gap_ext_ + std::max(Q_[i][j-1], R_[i][j-1] + gap_open_);
          if (Q_[i][j] == gap_ext_ + Q_[i][j-1]) hori_left_half_[i][j-1] = 1;
          if (Q_[i][j] == gap_ext_ + gap_open_ + R_[i][j-1]) hori_right_half_[i][j-1] = 1;
        }
        if (i != 0 && j != 0 ) {
          R_[i][j] = std::max(R_[i-1][j-1] + DiagScore_(ref_[j-1], query_[i-1]), std::max(Q_[i][j], P_[i][j]));
          if (R_[i][j]  == R_[i-1][j-1] + DiagScore_(ref_[j-1], query_[i-1])) diag_whole_[i][j] = 1;
        }
        if (R_[i][j]  == P_[i][j]) vert_whole_[i][j] = 1;
        if (R_[i][j]  == Q_[i][j]) hori_whole_[i][j] = 1;
      }
    }

    //after edge assignment bit array matrics
    for (int i = nrow_ - 1; i >= 0; --i) {
      for (int j = ncol_ - 1; j >= 0; --j) {
        if ((vert_whole_[i+1][j] == 0 || vert_bottom_half_[i][j] == 0) &&
            (hori_whole_[i][j+1] == 0 || hori_right_half_[i][j] == 0) &&
            diag_whole_[i+1][j+1] == 0) {
          vert_whole_[i][j] = 0;
          hori_whole_[i][j] = 0;
          diag_whole_[i][j] = 0;
        }
        if (vert_whole_[i+1][j] == 0  &&
            hori_whole_[i][j+1] == 0  &&
            diag_whole_[i+1][j+1] == 0) {
          continue;
        } else {
          if ( vert_whole_[i+1][j] == 1 && vert_top_half_[i][j] == 1) {
            vert_top_half_[i+1][j] = 1 -  vert_bottom_half_[i][j];
            vert_bottom_half_[i][j] = 1 - vert_whole_[i][j];
            vert_whole_[i][j] = 1;
          } else {
            vert_top_half_[i+1][j] = 0;
            vert_bottom_half_[i][j] = 0;
          }

          if ( hori_whole_[i][j + 1] == 1 && hori_left_half_[i][j] == 1) {
            hori_left_half_[i][j+1] = 1 -  hori_right_half_[i][j];
            hori_right_half_[i][j] = 1 - hori_whole_[i][j];
            hori_whole_[i][j] = 1;
          } else {
            hori_left_half_[i][j+1] = 0;
            hori_right_half_[i][j] = 0;
          }
        }
      }
    }
    // backtrack by bit array matrics
    int i = nrow_ - 1;
    int j = ncol_ - 1;
    bool must_go_vert = false;
    bool must_go_hori = false;
    while ( i != 0 || j != 0) {
      if (must_go_vert) {
        assert(vert_whole_[i][j] == 1);
        align_path_.GoNext(i-1, j);
        must_go_vert = vert_top_half_[i][j] == 1 ? 1 : 0;
        --i;
        continue;
      }
      if (must_go_hori) {
        assert(hori_whole_[i][j] == 1);
        align_path_.GoNext(i, j - 1);
        must_go_hori = hori_left_half_[i][j] == 1 ? 1 : 0;
        --j;
        continue;
      }

      if (diag_whole_[i][j] == 1) {
        align_path_.GoNext(i - 1 ,j - 1);
        --i;
        --j;
        //std::cerr << "diag i,j " << i <<", " << j <<std::endl;
      }
      else if (vert_whole_[i][j] == 1) {
        align_path_.GoNext(i-1, j);
        if (vert_top_half_[i][j] == 1) must_go_vert = true;
        --i;
        //std::cerr << "vert i,j " << i <<", " << j <<std::endl;
      } else if (hori_whole_[i][j] == 1) {
        align_path_.GoNext(i, j -1);
        if (hori_left_half_[i][j] == 1) must_go_hori = true;
        --j;
        //std::cerr << "hori i,j " << i <<", " << j <<std::endl;
      } else {
        assert (false);
      }
    }
  }

  void Print() const {
//    std::cerr << "R\n";
//    std::cerr << R_;
//    std::cerr << "P\n";
//    std::cerr << P_;
//    std::cerr << "Q\n";
//    std::cerr << Q_;
//    std::cerr << "a\n";
//    std::cerr << vert_whole_; // a
//    std::cerr << "b\n";
//    std::cerr << hori_whole_; // b
//    std::cerr << "c\n";
//    std::cerr << diag_whole_; // c
//    std::cerr << "d\n";
//    std::cerr << vert_top_half_; // d
//    std::cerr << "e\n";
//    std::cerr << vert_bottom_half_; // e
//    std::cerr << "f\n";
//    std::cerr << hori_left_half_; // f
//    std::cerr << "g\n";
//    std::cerr << hori_right_half_; // g
    std::string ref;
    std::string query;
    for(int i = align_path_.Ref().length() - 1; i >= 0; --i) {
      ref += align_path_.Ref()[i];
      query += align_path_.Alt()[i];
    }
    std::cout << ref << std::endl;
    std::cout << query << std::endl;
  }
};
int main() {
  AffineGap ag("AAAGGG", "TTAAAAGGGGTT");
  ag.Print();
  return 0;
}

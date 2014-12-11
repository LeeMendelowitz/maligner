/*
Test the align.cc methods
*/


#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

// Test the ChunkDatabase
#include "gtest/gtest.h"
#include "types.h"

#include "ScoreMatrix.h"
#include "align.h"
#include "utils.h"

PartialSums dummy_partial_sums;

// To use a test fixture, derive a class from testing::Test.
class AlignTest : public testing::Test {
 protected:  // You should make the members protected s.t. they can be
             // accessed from sub-classes.

  // virtual void SetUp() will be called before each test is run.  You
  // should define it if you need to initialize the varaibles.
  // Otherwise, this can be skipped.
  virtual void SetUp()
  {
    query_ = {1, 3, 5 ,9, 7, 9, 7, 5, 4, 1};
    ref_ = {4, 1, 3, 5 ,9, 7, 9, 7, 5, 4, 1, 9};
    size_t m = query_.size();
    size_t n = ref_.size();
    mat_ = ScoreMatrix(m+1, n+1);
  }

  // virtual void TearDown() will be called after each test is run.
  // You should define it if there is cleanup work to do.  Otherwise,
  // you don't have to provide it.
  //
  virtual void TearDown() {
  }

  // Declares the variables your tests want to use.
  ScoreMatrix mat_;
  IntVec query_;
  IntVec ref_;

};


TEST_F(AlignTest, align) {
    AlignOpts align_opts = AlignOpts(5.0, 3.0, 3, 3, 25);
    AlignTask align_task(query_, ref_, dummy_partial_sums, dummy_partial_sums, &mat_, align_opts);
    fill_score_matrix(align_task);
    cerr << "Done filling score matrix!";
    ASSERT_TRUE(true);
}


TEST_F(AlignTest, align_extra_rows) {

    // Align using a score matrix with extra rows.
    AlignOpts align_opts = AlignOpts(5.0, 3.0, 3, 3, 25);
    int m = query_.size() + 1;
    int n = ref_.size() + 1;

    ScoreMatrix * mat = new ScoreMatrix(m, n);
    ScoreMatrix * mat2 = new ScoreMatrix(m+100, n);

    AlignTask align_task(query_, ref_, dummy_partial_sums, dummy_partial_sums, mat, align_opts);
    //AlignTask align_task2(query_, ref_, mat2, align_opts);

    fill_score_matrix(align_task);
    //ill_score_matrix(align_task2);

    // Now build the best alignment
    ScoreCellPVec trail, trail2;

    bool result = get_best_alignment(align_task, trail);
    //bool result2 = get_best_alignment(align_task2, trail2);

    std::cout << "Trail1: \n" << trail << "\n";
    //std::cout << "Trail2: \n" << trail2 << "\n";

    ASSERT_TRUE(result);
    //ASSERT_TRUE(result2);

    ChunkVec query_chunks, query_chunks2, ref_chunks, ref_chunks2;

    Alignment * a = alignment_from_trail(align_task, trail);
    //Alignment a2 = alignment_from_trail(align_task2, trail2);

    // Open output file for outputing score matrix
    ofstream fout("score_matrix.txt");
    fout << *mat << "\n";
    fout.close();

    delete a;

}


TEST_F(AlignTest, make_align) {

    // Align using a score matrix with extra rows.
    AlignOpts align_opts = AlignOpts(5.0, 3.0, 3, 3, 25);
    IntVec query = {1, 3, 5 ,9, 17, 7, 5, 4, 1};
    IntVec ref = {4, 1, 3, 5 ,9, 7, 9, 7, 5, 4, 1, 9};

    int m = query.size() + 1;
    int n = ref.size() + 1;

    ScoreMatrix mat = ScoreMatrix(m, n);
    AlignTask align_task(query, ref, dummy_partial_sums, dummy_partial_sums, &mat, align_opts);

    fill_score_matrix(align_task);

    // Now build the best alignment
    ScoreCellPVec trail;
    bool result = get_best_alignment(align_task, trail);
    ASSERT_TRUE(result);
    
    Alignment * a(alignment_from_trail(align_task, trail));
    std::cerr << "Made alignment: " << *a;
    delete a;
}

TEST_F(AlignTest, partial_sum_test) {
  IntVec q = {1,2,3,4,5};
  const int missed_sites = 2;
  vector<IntVec> partial_sum = make_partial_sums(q, missed_sites);

  ASSERT_TRUE(partial_sum.size() == q.size());

  for (size_t i = 0; i < partial_sum.size(); i++) {
    ASSERT_TRUE(partial_sum[i].size() == missed_sites+1);
    cout << "i: " << i << ", partial_sum: " << partial_sum[i] << "\n";
  }
}


TEST_F(AlignTest, partial_sum_align_test) {
  
  // Align using a score matrix with extra rows.
  const int missed_sites = 3;
  AlignOpts align_opts = AlignOpts(5.0, 3.0, missed_sites, missed_sites, 25);

  IntVec query = {1, 3, 5 ,9, 17, 7, 5, 4, 1};
  IntVec ref = {4, 1, 3, 5 ,9, 7, 9, 7, 5, 4, 1, 9};

  PartialSums query_partial_sum = make_partial_sums(query, missed_sites);
  PartialSums ref_partial_sum = make_partial_sums(ref, missed_sites);

  int m = query.size() + 1;
  int n = ref.size() + 1;

  ScoreMatrix mat = ScoreMatrix(m, n);
  AlignTask align_task(query, ref, query_partial_sum, ref_partial_sum, &mat, align_opts);

  fill_score_matrix_using_partials(align_task);

  // Now build the best alignment
  ScoreCellPVec trail;
  bool result = get_best_alignment(align_task, trail);
  ASSERT_TRUE(result);
  
  Alignment * a(alignment_from_trail(align_task, trail));
  std::cerr << "Made alignment: " << *a;
  delete a;
}

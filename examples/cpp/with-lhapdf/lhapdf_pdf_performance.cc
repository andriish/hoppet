#include "hoppet.h"
#include "LHAPDF/LHAPDF.h"
#include "lhapdf_c_api.h"

#include <chrono>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

using namespace std;

static void benchmark_pdf_performance(const string& pdfname, int member = 0) {
  unique_ptr<LHAPDF::PDF> lhapdf_pdf(LHAPDF::mkPDF(pdfname, member));
  unique_ptr<HOPPETPDF::PDF> hoppet_pdf(HOPPETPDF::mkPDF(pdfname, member));

  if (!lhapdf_pdf || !hoppet_pdf) {
    cerr << "Failed to load PDF set: " << pdfname << endl;
    return;
  }

  const int npoints = 800;
  vector<double> xvals(npoints), qvals(npoints);
  const double xmin = 1e-5;
  const double xmax = 0.9;
  const double qmin = 1.5;
  const double qmax = 7000.0;

  for (int i = 0; i < npoints; ++i) {
    xvals[i] = exp(log(xmin) + i * (log(xmax) - log(xmin)) / (npoints - 1));
    qvals[i] = exp(log(qmin) + i * (log(qmax) - log(qmin)) / (npoints - 1));
  }

  vector<double> lhapdf_vec(13, 0.0);
  vector<double> hoppet_vec(13, 0.0);

  // Warm-up evaluation to initialize caches / reduce first-call overhead.
  lhapdf_pdf->xfxQ2(xvals.front(), qvals.front() * qvals.front(), lhapdf_vec);
  hoppet_pdf->xfxQ2(xvals.front(), qvals.front() * qvals.front(), hoppet_vec);

  auto t1 = chrono::high_resolution_clock::now();
  for (double x : xvals) {
    for (double q : qvals) {
      lhapdf_pdf->xfxQ2(x, q * q, lhapdf_vec);
    }
  }
  auto t2 = chrono::high_resolution_clock::now();
  const double lh_time_ns = chrono::duration<double, nano>(t2 - t1).count() / (npoints * npoints);

  t1 = chrono::high_resolution_clock::now();
  for (double x : xvals) {
    for (double q : qvals) {
      hoppet_pdf->xfxQ2(x, q * q, hoppet_vec);
    }
  }
  t2 = chrono::high_resolution_clock::now();
  const double hp_time_ns = chrono::duration<double, nano>(t2 - t1).count() / (npoints * npoints);

  cout << "LHAPDF::PDF::xfxQ2: " << lh_time_ns << " ns/eval" << endl;
  cout << "HOPPETPDF::PDF::xfxQ2: " << hp_time_ns << " ns/eval" << endl;
  cout << "Speedup factor (LHAPDF/HOPPET): " << (lh_time_ns / hp_time_ns) << "x" << endl;
}

int main(int argc, char* argv[]) {
  string pdfname = "PDF4LHC21_40";
  if (argc > 1) {
    pdfname = argv[1];
  }

  cout << "Comparing xfxQ2 performance for LHAPDF::PDF vs HOPPETPDF::PDF using set: " << pdfname << endl;
  benchmark_pdf_performance(pdfname, 0);
  return 0;
}

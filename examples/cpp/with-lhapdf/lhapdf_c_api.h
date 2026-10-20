#ifndef LHAPDF_C_API_H
#define LHAPDF_C_API_H

#include <cstdint>
#ifdef __cplusplus
#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#endif

#ifdef __cplusplus
extern "C" {
#endif

struct pdf_grid_t {
  int n;
  int nQ;
  int nFlav;
  double* x_ptr;
  double* q_ptr;
  double* f_ptr;
  int* flav_ptr;
};

struct pdf_location_t {
  char pdfsets_path[1024];
  char pdfdata_path[1024];
  char pdfinfo_path[1024];
  int valid;
};

struct pdf_metadata_t {
  char SetDesc[1024];
  int SetIndex;
  char Authors[1024];
  char Reference[1024];
  char Format[1024];
  char DataVersion[1024];
  int NumMembers;
  int Particle;
  intptr_t Flavors;
  int OrderQCD;
  char FlavorScheme[1024];
  int NumFlavors;
  char ErrorType[1024];
  double XMin;
  double XMax;
  double QMin;
  double QMax;
  double MZ;
  double MUp;
  double MDown;
  double MStrange;
  double MCharm;
  double MBottom;
  double MTop;
  double AlphaS_MZ;
  int AlphaS_OrderQCD;
  char AlphaS_Type[1024];
  intptr_t AlphaS_Qs;
  intptr_t AlphaS_Vals;
};

pdf_grid_t pdf_read_grid_c(const char* fname, int flen);
pdf_location_t pdf_locate_c(const char* pdfname, int pdfname_len);
pdf_metadata_t pdf_read_metadata_c(const char* fname, int flen);
void pdf_get_values_c(pdf_grid_t data, double x, double Q, int flv, double* val);
double pdf_get_alphas_c(pdf_metadata_t meta, double Q);

#ifdef __cplusplus
}

namespace HOPPETPDF {
class InitialPDF {
public:
  InitialPDF(const std::string& pdfname, int member = 0)
      : loc_(),
        member_(member),
        metadata_(),
        grid_() {
    loc_ = pdf_locate_c(pdfname.c_str(), static_cast<int>(pdfname.size()));
    std::cout << "[DEBUG] HOPPETPDF::InitialPDF::locate_pdf returned info_path='" << loc_.pdfinfo_path << "'" << std::endl;
    std::cout << "[DEBUG] HOPPETPDF::InitialPDF::locate_pdf returned dat_path='" << loc_.pdfdata_path << "'" << std::endl;

    std::cout << "[DEBUG] Reading PDF info from: " << loc_.pdfinfo_path << std::endl;
    metadata_ = pdf_read_metadata_c(loc_.pdfinfo_path, static_cast<int>(std::strlen(loc_.pdfinfo_path)));

    std::string data_path = loc_.pdfdata_path;
    if (member_ != 0) {
      const std::size_t pos = data_path.rfind("_0000.dat");
      if (pos != std::string::npos) {
        std::ostringstream oss;
        oss << '_' << std::setw(4) << std::setfill('0') << member_ << ".dat";
        data_path.replace(pos, 9, oss.str());
      }
    }

    std::cout << "[DEBUG] Reading PDF data from: " << data_path << std::endl;
    grid_ = pdf_read_grid_c(data_path.c_str(), static_cast<int>(data_path.size()));
    std::cout << "[DEBUG] HOPPETPDF::InitialPDF::read_grid returned n=" << grid_.n << " nQ=" << grid_.nQ
              << " nFlav=" << grid_.nFlav << std::endl;
  }

  ~InitialPDF() {
    std::memset(&loc_, 0, sizeof(loc_));
    std::memset(&metadata_, 0, sizeof(metadata_));
    grid_.n = 0;
    grid_.nQ = 0;
    grid_.nFlav = 0;
    grid_.x_ptr = nullptr;
    grid_.q_ptr = nullptr;
    grid_.f_ptr = nullptr;
    grid_.flav_ptr = nullptr;
  }

  int orderQCD() const { return metadata_.OrderQCD; }
  double xMin() const { return metadata_.XMin; }
  double xMax() const { return metadata_.XMax; }
  double qMin() const { return metadata_.QMin; }
  double qMax() const { return metadata_.QMax; }

  double quarkMass(int flavor) const {
    switch (flavor) {
      case 1: return metadata_.MUp;
      case 2: return metadata_.MDown;
      case 3: return metadata_.MStrange;
      case 4: return metadata_.MCharm;
      case 5: return metadata_.MBottom;
      case 6: return metadata_.MTop;
      default: return 0.0;
    }
  }

  double alphasQ(double q) const {
    return pdf_get_alphas_c(metadata_, q);
  }

  double xfxQ2(int id, double x, double q2) const {
    double val = 0.0;
    pdf_get_values_c(grid_, x, std::sqrt(q2), id, &val);
    return val;
  }

  void xfxQ2(double x, double q2, std::vector<double>& out) const {
    out.assign(13, 0.0);
    for (int flv = -6; flv <= 6; ++flv) {
      out[flv + 6] = xfxQ2(flv, x, q2);
    }
  }

  bool hasFlavor(int flavor) const {
    if (grid_.nFlav <= 0 || grid_.flav_ptr == nullptr) {
      return false;
    }
    const int* flav = grid_.flav_ptr;
    for (int i = 0; i < grid_.nFlav; ++i) {
      if (flav[i] == flavor) {
        return true;
      }
    }
    return false;
  }

  const pdf_grid_t& grid() const { return grid_; }
  const pdf_metadata_t& metadata() const { return metadata_; }

private:
  pdf_location_t loc_;
  pdf_metadata_t metadata_;
  int member_;
  pdf_grid_t grid_;
};

class PDF {


  private:
/// Interface that supplies PDF values from the wrapper-backed .dat grid.
/// The flavour ordering matches LHAPDF's convention (-6 … +6).
inline static std::shared_ptr<InitialPDF> g_initial_pdf = nullptr;
static void initial_interface(const double & x, const double & Q, double * res) {
  std::vector<double> vals(13, 0.0);
  g_initial_pdf->xfxQ2(x, Q * Q, vals);
  std::copy(vals.begin(), vals.end(), res);
}
  public:
  PDF(const std::string& pdfname, int member = 0){
  auto initial_pdf = std::make_shared<InitialPDF>(pdfname, member);
  metadata_ = initial_pdf->metadata();
  member_ = member;

  const double xmin = metadata_.XMin;
  const double xmax = metadata_.XMax;
  const double Qmin = metadata_.QMin;
  const double Qmax = metadata_.QMax;
  const int nloop = metadata_.OrderQCD + 1;
  const double mc = metadata_.MCharm;
  const double mb = metadata_.MBottom;
  double mt = metadata_.MTop;
  if (mt <= 0.0) {
    mt = 2.0 * Qmax;
  }

  const double ymax = static_cast<double>(std::ceil(std::log(1.0 / xmin)));
  double dy = 0.05;
  double dlnlnQ = dy / 4.0;
  if (ymax > 15.0) {
    dlnlnQ = dy / 8.0;
  }
  const int order = -6;
  const int yorder = 2;
  const int lnlnQorder = 2;

  hoppetSetPoleMassVFN(mc, mb, mt);
  hoppetSetYLnlnQInterpOrders(yorder, lnlnQorder);
  hoppetStartExtended(ymax, dy, Qmin, Qmax, dlnlnQ, nloop, order, hoppet::factscheme_MSbar);
  const double Q0 = Qmin;
  hoppetSetCoupling(initial_pdf->alphasQ(Q0), Q0, nloop);
  g_initial_pdf = initial_pdf;
  hoppetAssign(initial_interface);
  g_initial_pdf = nullptr;
  }
  
  int orderQCD() const { return metadata_.OrderQCD; }
  double xMin() const { return metadata_.XMin; }
  double xMax() const { return metadata_.XMax; }
  double qMin() const { return metadata_.QMin; }
  double qMax() const { return metadata_.QMax; }

  double quarkMass(int flavor) const {
    switch (flavor) {
      case 1: return metadata_.MUp;
      case 2: return metadata_.MDown;
      case 3: return metadata_.MStrange;
      case 4: return metadata_.MCharm;
      case 5: return metadata_.MBottom;
      case 6: return metadata_.MTop;
      default: return 0.0;
    }
  }

  double alphasQ(double q) const {
    return hoppetAlphaS(q);
  }

  double xfxQ2(int id, double x, double q2) const {
    if (id < -6 || id > 6) {
      return 0.0;
    }
    double hoppetpdf[13];
    // Equivalent hoppet call
    hoppetEval(x, std::sqrt(q2), hoppetpdf);
    return hoppetpdf[id + 6];
  }

  void xfxQ2(double x, double q2, std::vector<double>& out) const {
    out.assign(13, 0.0);
    double hoppetpdf[13];
    // Equivalent hoppet call
    hoppetEval(x, std::sqrt(q2), hoppetpdf);
    out.assign(hoppetpdf, hoppetpdf + 13);
  }
  
  ~PDF() {
  }
private:
  pdf_metadata_t metadata_;
  int member_;

};
PDF* mkPDF(const std::string& pdfname, int member = 0) {
  return new PDF(pdfname, member);
}
InitialPDF* mkInitialPDF(const std::string& pdfname, int member = 0) {
  return new InitialPDF(pdfname, member);
}


} // namespace HOPPETPDF

#endif

#endif

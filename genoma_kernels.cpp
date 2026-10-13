#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <unordered_map>

namespace py = pybind11;

inline char complement_char(char c) {
    switch (c) {
        case 'A': case 'a': return 'T';
        case 'C': case 'c': return 'G';
        case 'G': case 'g': return 'C';
        case 'T': case 't': return 'A';
        case 'N': case 'n': return 'N';
        default: return c;
    }
}

std::string reverse_complement_cpp(const std::string& seq) {
    std::string rev(seq.rbegin(), seq.rend());
    for (char& c : rev) c = complement_char(c);
    return rev;
}

double puntuacion_pwm_cpp(const std::string& ventana, 
                          const std::vector<std::unordered_map<char, double>>& matrix, 
                          double min_score, double max_score) {
    if (ventana.size() != matrix.size()) return -999.0;
    double score = 0.0;
    for (size_t i = 0; i < ventana.size(); ++i) {
        char b = std::toupper(ventana[i]);
        auto it = matrix[i].find(b);
        if (it == matrix[i].end()) return -999.0;
        score += it->second;
    }
    if (max_score <= min_score) return 0.0;
    return (score - min_score) / (max_score - min_score);
}

py::dict mejor_par_pwm_centrada_cpp(const std::string& seq_anc, 
                                     const std::string& seq_hum, 
                                     int centro, 
                                     const std::vector<std::unordered_map<char, double>>& matrix,
                                     double min_score, double max_score, int width) {
    int lo = std::max(0, centro - width + 1);
    int hi = std::min(centro, (int)seq_anc.length() - width);
    
    bool found = false;
    double best_diff = -1.0;
    double best_max_score = -1.0;
    double best_sa = 0.0, best_sh = 0.0;
    int best_start = 0;
    std::string best_orient = "+";
    std::string best_aa = "", best_hh = "";

    for (int start = lo; start <= hi; ++start) {
        std::string va = seq_anc.substr(start, width);
        std::string vh = seq_hum.substr(start, width);

        for (const std::string& orient : {"+", "-"}) {
            std::string aa = va, hh = vh;
            if (orient == "-") {
                aa = reverse_complement_cpp(va);
                hh = reverse_complement_cpp(vh);
            }

            double sa = puntuacion_pwm_cpp(aa, matrix, min_score, max_score);
            double sh = puntuacion_pwm_cpp(hh, matrix, min_score, max_score);

            if (sa < -100.0 || sh < -100.0) continue;

            double diff = std::abs(sa - sh);
            double m_score = std::max(sa, sh);

            if (!found || diff > best_diff || (diff == best_diff && m_score > best_max_score)) {
                found = true;
                best_diff = diff;
                best_max_score = m_score;
                best_sa = sa;
                best_sh = sh;
                best_start = start;
                best_orient = orient;
                best_aa = aa;
                best_hh = hh;
            }
        }
    }

    if (!found) return py::dict();

    py::dict result;
    result["score_ancestral"] = best_sa;
    result["score_humano"] = best_sh;
    result["start"] = best_start;
    result["orientacion"] = best_orient;
    result["secuencia_ancestral"] = best_aa;
    result["secuencia_humana"] = best_hh;
    return result;
}

PYBIND11_MODULE(genoma_kernels, m) {
    m.def("puntuacion_pwm_cpp", &puntuacion_pwm_cpp);
    m.def("mejor_par_pwm_centrada_cpp", &mejor_par_pwm_centrada_cpp);
    m.def("reverse_complement_cpp", &reverse_complement_cpp);
}
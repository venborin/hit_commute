#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

using Eigen::MatrixXd;
using Eigen::VectorXd;

struct Residue {
    std::string resName;
    char chain;        // ' ' if missing
    int resSeq;        // residue sequence number
    char iCode;        // insertion code or ' '
    Eigen::Vector3d ca; // CA coordinates
    int index;         // 0..N-1
};

struct Args {
    std::string pdbPath;
    std::optional<std::string> targetResidue; // e.g., "A:42" or "42" or "A:42B"
    double cutoff = 7.5;
    bool matrixMode = false; // print full matrices if no --res
};

static void die(const std::string& msg) {
    std::cerr << "Error: " << msg << "\n";
    std::exit(1);
}

// Trim helpers
static std::string trim(const std::string& s) {
    size_t b = s.find_first_not_of(" \t\r\n");
    if (b == std::string::npos) return "";
    size_t e = s.find_last_not_of(" \t\r\n");
    return s.substr(b, e - b + 1);
}

static bool starts_with(const std::string& s, const std::string& p) {
    return s.size() >= p.size() && std::equal(p.begin(), p.end(), s.begin());
}

// Parse residue selector like "A:42", "A:42B", or "42"
struct ResidueSelector {
    std::optional<char> chain; // absent => no chain constraint
    int resSeq = 0;
    char iCode = ' '; // optional insertion code
};
static ResidueSelector parseResidueSelector(const std::string& sel) {
    ResidueSelector rs;
    std::string s = trim(sel);
    auto colon = s.find(':');
    if (colon == std::string::npos) {
        // just a number (no chain)
        rs.chain = std::nullopt;
        // check trailing letter as insertion code
        // e.g., "42B"
        size_t pos = s.find_first_not_of("0123456789");
        if (pos == std::string::npos) {
            rs.resSeq = std::stoi(s);
            rs.iCode = ' ';
        } else {
            rs.resSeq = std::stoi(s.substr(0, pos));
            if (pos == s.size()-1 && std::isalpha(static_cast<unsigned char>(s.back())))
                rs.iCode = s.back();
            else
                die("Bad residue selector: " + sel);
        }
    } else {
        if (colon == 0) die("Bad residue selector (missing chain): " + sel);
        if (colon > 1) die("Chain must be a single character before ':' in selector: " + sel);
        rs.chain = s[0];
        std::string rest = s.substr(colon+1);
        // rest can be "42" or "42B"
        size_t pos = rest.find_first_not_of("0123456789");
        if (pos == std::string::npos) {
            rs.resSeq = std::stoi(rest);
            rs.iCode = ' ';
        } else {
            rs.resSeq = std::stoi(rest.substr(0, pos));
            if (pos == rest.size()-1 && std::isalpha(static_cast<unsigned char>(rest.back())))
                rs.iCode = rest.back();
            else
                die("Bad residue selector: " + sel);
        }
    }
    return rs;
}

static void usage() {
    std::cout <<
R"(protein_ht — hitting & commute times on residue interaction networks

Usage:
  protein_ht --pdb FILE.pdb [--cutoff Å] [--res SEL]

Options:
  --pdb PATH         Input PDB file.
  --cutoff FLOAT     Cα–Cα cutoff in Å for edges (default 7.5).
  --res SEL          Compute times from/to this residue only.
                     SEL formats: "A:42", "A:42B", or "42" (no chain).
Notes:
  - If the PDB lacks chain IDs, residues have blank chain (' '), so use --res 42.
  - Output uses the corrected formulas for hitting/commute times via L^+.
)";
}

static Args parseArgs(int argc, char** argv) {
    Args a;
    if (argc < 3) { usage(); die("Not enough arguments."); }
    for (int i=1; i<argc; ++i) {
        std::string t = argv[i];
        if (t == "--pdb") {
            if (++i >= argc) die("Missing value for --pdb");
            a.pdbPath = argv[i];
        } else if (t == "--res") {
            if (++i >= argc) die("Missing value for --res");
            a.targetResidue = argv[i];
        } else if (t == "--cutoff") {
            if (++i >= argc) die("Missing value for --cutoff");
            a.cutoff = std::stod(argv[i]);
        } else if (t == "--help" || t == "-h") {
            usage(); std::exit(0);
        } else {
            die("Unknown argument: " + t);
        }
    }
    if (a.pdbPath.empty()) die("Please provide --pdb PATH");
    a.matrixMode = !a.targetResidue.has_value();
    return a;
}

// Parse ATOM records for CA atoms; tolerate missing chain (blank)
static std::vector<Residue> readResiduesFromPDB(const std::string& path) {
    std::ifstream in(path);
    if (!in) die("Cannot open PDB: " + path);

    // Key by (chain, resSeq, iCode)
    struct Key { char chain; int resseq; char icode;
        bool operator<(const Key& o) const {
            return std::tie(chain, resseq, icode) < std::tie(o.chain, o.resseq, o.icode);
        }
    };

    std::map<Key, Residue> residues;
    std::string line;

    while (std::getline(in, line)) {
        if (line.size() < 54) continue;
        if (!(starts_with(line, "ATOM  ") || starts_with(line, "HETATM"))) continue;

        // Columns per PDB v3.3+:
        //  1-6  record name
        // 13-16 atom name
        // 17    altLoc
        // 18-20 resName
        // 22    chainID
        // 23-26 resSeq
        // 27    iCode
        // 31-38 x, 39-46 y, 47-54 z  (1-based)
        std::string atomName = (line.size() >= 16) ? line.substr(12,4) : "";
        atomName = trim(atomName);
        if (atomName != "CA") continue;

        char altLoc = (line.size() >= 17) ? line[16] : ' ';
        if (altLoc != ' ' && altLoc != 'A') continue; // ignore alt confs except blank/A

        std::string resName = (line.size() >= 20) ? trim(line.substr(17,3)) : "UNK";
        char chain = (line.size() >= 22) ? line[21] : ' ';
        if (chain == '\0') chain = ' '; // normalize

        int resSeq = 0;
        try {
            resSeq = std::stoi((line.size() >= 26) ? line.substr(22,4) : "0");
        } catch (...) {
            continue;
        }
        char iCode = (line.size() >= 27) ? line[26] : ' ';
        if (iCode == '\0') iCode = ' ';

        double x = std::stod(line.substr(30,8));
        double y = std::stod(line.substr(38,8));
        double z = std::stod(line.substr(46,8));

        Key k{chain, resSeq, iCode};
        // If multiple CA entries exist (rare), keep the first.
        if (residues.find(k) == residues.end()) {
            Residue r;
            r.resName = resName;
            r.chain = chain;
            r.resSeq = resSeq;
            r.iCode = iCode;
            r.ca << x, y, z;
            r.index = -1;
            residues[k] = r;
        }
    }

    std::vector<Residue> out;
    out.reserve(residues.size());
    int idx = 0;
    for (auto &kv : residues) {
        kv.second.index = idx++;
        out.push_back(kv.second);
    }
    if (out.empty()) die("No residues (Cα) parsed from PDB.");
    return out;
}

// Build adjacency (unweighted) with cutoff
static MatrixXd buildAdjacency(const std::vector<Residue>& R, double cutoff) {
    const int n = (int)R.size();
    MatrixXd A = MatrixXd::Zero(n, n);
    const double c2 = cutoff * cutoff;
    for (int i=0; i<n; ++i) {
        for (int j=i+1; j<n; ++j) {
            double d2 = (R[i].ca - R[j].ca).squaredNorm();
            if (d2 <= c2) {
                A(i,j) = 1.0;
                A(j,i) = 1.0;
            }
        }
    }
    return A;
}

// Laplacian L = D - A
static MatrixXd buildLaplacian(const MatrixXd& A, VectorXd& degree, double& volume) {
    const int n = (int)A.rows();
    degree = A.rowwise().sum();
    volume = degree.sum(); // sum of degrees (2 * number of edges)
    MatrixXd L = -A;
    for (int i=0; i<n; ++i) L(i,i) = degree(i);
    return L;
}

// Pseudoinverse of symmetric PSD Laplacian via eigen-decomposition
static MatrixXd laplacianPseudoinverse(const MatrixXd& L) {
    Eigen::SelfAdjointEigenSolver<MatrixXd> es(L);
    if (es.info() != Eigen::Success) die("Eigen decomposition failed.");

    const double eps = 1e-10; // threshold for zero eigenvalue(s)
    VectorXd evals = es.eigenvalues();
    MatrixXd U = es.eigenvectors();

    MatrixXd invDiag = MatrixXd::Zero(L.rows(), L.cols());
    for (int i=0; i<evals.size(); ++i) {
        double lam = evals(i);
        if (lam > eps) invDiag(i,i) = 1.0 / lam;
        else invDiag(i,i) = 0.0; // project out nullspace (the all-ones vector)
    }
    MatrixXd Lplus = U * invDiag * U.transpose();
    return Lplus;
}

static bool isConnectedFromSpectrum(const MatrixXd& L) {
    Eigen::SelfAdjointEigenSolver<MatrixXd> es(L);
    if (es.info() != Eigen::Success) return false;
    const double eps = 1e-10;
    int zeros = 0;
    for (int i=0; i<es.eigenvalues().size(); ++i) {
        if (es.eigenvalues()(i) < eps) ++zeros;
    }
    return zeros == 1;
}

static std::string residLabel(const Residue& r) {
    // show '-' for blank chain to make it clear in output
    char chainOut = (r.chain == ' ' ? '-' : r.chain);
    std::ostringstream oss;
    oss << chainOut << ":" << r.resSeq;
    if (r.iCode != ' ') oss << r.iCode;
    oss << " (" << r.resName << ")";
    return oss.str();
}

static std::optional<int> findResidueIndex(
    const std::vector<Residue>& R, const ResidueSelector& sel)
{
    for (const auto& r : R) {
        bool chainMatch = true;
        if (sel.chain.has_value()) {
            chainMatch = (r.chain == sel.chain.value());
        }
        if (chainMatch && r.resSeq == sel.resSeq && r.iCode == sel.iCode) {
            return r.index;
        }
    }
    return std::nullopt;
}

int main(int argc, char** argv) {
    Args args = parseArgs(argc, argv);

    auto residues = readResiduesFromPDB(args.pdbPath);
    const int n = (int)residues.size();

    auto A = buildAdjacency(residues, args.cutoff);
    VectorXd degree;
    double volume = 0.0;
    auto L = buildLaplacian(A, degree, volume);

    if (!isConnectedFromSpectrum(L)) {
        die("Graph is disconnected. Increase --cutoff or restrict to a connected chain.");
    }

    auto Lplus = laplacianPseudoinverse(L);

    // Hitting/Commute formulas (equivalent to corrected expressions):
    // H_{i->j} = vol * (L^+_{jj} - L^+_{ij})
    // C_{ij}   = vol * (L^+_{ii} + L^+_{jj} - 2 L^+_{ij})

    if (args.matrixMode) {
        std::cout << "# N=" << n << ", cutoff=" << args.cutoff
                  << ", vol=" << volume << "\n";
        std::cout << "# Residue index map:\n";
        for (const auto& r : residues) {
            std::cout << r.index << "  " << residLabel(r) << "\n";
        }

        std::cout << "\n# Commute time matrix (N x N)\n";
        for (int i=0; i<n; ++i) {
            for (int j=0; j<n; ++j) {
                double cij = volume * (Lplus(i,i) + Lplus(j,j) - 2.0*Lplus(i,j));
                std::cout << (j? " " : "") << cij;
            }
            std::cout << "\n";
        }

        std::cout << "\n# Hitting time matrix H_{i->j} (rows: source i, cols: target j)\n";
        for (int i=0; i<n; ++i) {
            for (int j=0; j<n; ++j) {
                double hij = volume * (Lplus(j,j) - Lplus(i,j));
                std::cout << (j? " " : "") << hij;
            }
            std::cout << "\n";
        }
        return 0;
    }

    // Single-residue mode
    ResidueSelector sel = parseResidueSelector(*args.targetResidue);
    auto idxOpt = findResidueIndex(residues, sel);
    if (!idxOpt.has_value()) {
        std::ostringstream oss;
        oss << "Residue not found: " << *args.targetResidue
            << " (remember: blank chain residues use --res RESSEQ with no chain).";
        die(oss.str());
    }
    int s = idxOpt.value();

    std::cout << "# Target residue: " << residLabel(residues[s])
              << " (index " << s << ")\n";
    std::cout << "# N=" << n << ", cutoff=" << args.cutoff
              << ", vol=" << volume << "\n";
    std::cout << "index\tresidue\tH(s->j)\tH(j->s)\tCommute\n";

    for (int j=0; j<n; ++j) {
        double H_sj = volume * (Lplus(j,j) - Lplus(s,j));
        double H_js = volume * (Lplus(s,s) - Lplus(j,s));
        double C_sj = H_sj + H_js; // or volume*(L^+_{ss}+L^+_{jj}-2L^+_{sj})
        std::cout << j << "\t" << residLabel(residues[j]) << "\t"
                  << H_sj << "\t" << H_js << "\t" << C_sj << "\n";
    }

    return 0;
}

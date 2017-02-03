#include <jarngreipr/io/PDBReader.hpp>
#include <jarngreipr/model/ClementiGo.hpp>
#include <mjolnir/core/DefaultTraits.hpp>
#include <mjolnir/core/RandomNumberGenerator.hpp>
#include <map>

typedef mjolnir::DefaultTraits traits;

constexpr static double T = 300.0;
constexpr static double kB = 1.986231313e-3;

int main(int argc, char **argv)
{
    if(argc != 2)
    {
        std::cerr << "./jarngreipr [file.pdb]" << std::endl;
        return 1;
    }
    mjolnir::RandomNumberGenerator<traits> rng(1234567);
    std::map<std::string, double> mass;
    mass["ALA"] = 71.09;
    mass["ARG"] = 156.19;
    mass["ASN"] = 114.11;
    mass["ASP"] = 115.09;
    mass["CYS"] = 103.15;
    mass["GLN"] = 128.14;
    mass["GLU"] = 129.12;
    mass["GLY"] = 57.05;
    mass["HIS"] = 137.14;
    mass["ILE"] = 113.16;
    mass["LEU"] = 113.16;
    mass["LYS"] = 128.17;
    mass["MET"] = 131.19;
    mass["PHE"] = 147.18;
    mass["PRO"] = 97.12;
    mass["SER"] = 87.08;
    mass["THR"] = 101.11;
    mass["TRP"] = 186.21;
    mass["TYR"] = 163.18;
    mass["VAL"] = 99.14;
    std::cerr << "parameter setup end" << std::endl;

    // hardly-coding c-alpha model
    std::string fname(argv[1]);
    jarngreipr::PDBReader<traits> pdb;
    auto atoms  = pdb.read(fname);
    std::cerr << atoms.size() << " particles found" << std::endl;
    auto chains = pdb.parse(atoms);
    std::cerr << chains.size() << " chains found" << std::endl;

    jarngreipr::ClementiGo<traits> clementigo;
    clementigo.make(chains);

    // output

    // particles

    std::vector<double> fric_consts;
    std::cout << "# particles {{{" << std::endl;
    std::cout << "particles = [" << std::endl;
    for(auto iter = clementigo.beads().cbegin();
            iter != clementigo.beads().cend(); ++iter) // bead
    {
        const std::string name = (*iter)->name();
        const typename traits::real_type m = mass[name];

        const auto pos = (*iter)->position(0);
        const auto vel = rng.maxwell_boltzmann(m, T, kB);
        const auto acc = traits::coordinate_type(0., 0., 0.);

        std::cout << std::scientific;
        std::cout << "{mass = " << m << ", position = ["
                  << pos[0] << ", " << pos[1] << ", " << pos[2]
                  << "], velocity = ["
                  << vel[0] << ", " << vel[1] << ", " << vel[2]
                  << "]},"
                  << std::endl;

        const typename traits::real_type f = 0.005 * 168.7 / m;
        fric_consts.push_back(f);
    }
    std::cout << "]" << std::endl;
    std::cout << "# }}}" << std::endl;

    std::cout << "# friction {{{" << std::endl;
    std::cout << "friction_constant = [" << std::endl;
    for(auto iter = fric_consts.cbegin(); iter != fric_consts.cend(); ++iter)
        std::cout << *iter << ", " << std::endl;
    std::cout << "]" << std::endl;
    std::cout << "# }}}" << std::endl;

    clementigo.output_bond(std::cout);
    clementigo.output_go(std::cout);
    clementigo.output_angle(std::cout);
    clementigo.output_dihedral(std::cout);

    // exvs
    std::cout << "[[globalforcefield]]" << std::endl;
    std::cout << "interaction = \"Global\"" << std::endl;
    std::cout << "potential = \"ExcludedVolume\"" << std::endl;
    std::cout << "epsilon = 0.6" << std::endl;
    std::cout << "# params{{{" << std::endl;
    std::cout << "parameters = [" << std::endl;
    for(std::size_t i=0; i<clementigo.beads().size(); ++i)
    {
        std::cout << "{sigma = 2.0}," << std::endl;
    }
    std::cout << "]" << std::endl;
    std::cout << "# }}}" << std::endl;

    clementigo.output_exception(std::cout);

    return 0;
}

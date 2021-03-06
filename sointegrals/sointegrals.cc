#include <libplugin/plugin.h>
#include "psi4-dec.h"
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libmints/sointegral_twobody.h>
#include <libpsio/psio.h>

INIT_PLUGIN
using namespace boost;
namespace psi{ namespace sointegrals{
extern "C" int
read_options(std::string name, Options &options)
{
    if (name == "SOINTEGRALS"|| options.read_globals()) {
        /*- The amount of information printed
            to the output file -*/
        options.add_int("PRINT", 1);
        /*- Whether to compute two-electron integrals -*/
        options.add_bool("DO_TEI", true);
    }
    return true;
}
class ERIPrinter
{
public:
    // Our functor...the nice thing about using a C++ functor is that the
    // code here is inlined by the compiler.
    void operator() (int pabs, int qabs, int rabs, int sabs,
                     int pirrep, int pso,
                     int qirrep, int qso,
                     int rirrep, int rso,
                     int sirrep, int sso,
                     double value)
    {
        fprintf(outfile, "\t(%2d %2d | %2d %2d) = %20.10lf\n",
                        pabs, qabs, rabs, sabs, value);
    }
};
extern "C" PsiReturnType
sointegrals(Options &options)
{
    int print = options.get_int("PRINT");
    int doTei = options.get_bool("DO_TEI");
    shared_ptr<Molecule> molecule = Process::environment.molecule();
    // Form basis object:
    // Create a basis set parser object.
    shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    // Construct a new basis set.
    shared_ptr<BasisSet> aoBasis = BasisSet::construct(parser, molecule, "BASIS");
    // The integral factory oversees the creation of integral objects
    shared_ptr<IntegralFactory> integral(new IntegralFactory
            (aoBasis, aoBasis, aoBasis, aoBasis));
    // N.B. This should be called after the basis has been built, because
    // the geometry has not been
    // fully initialized until this time.
    molecule->print();
    // The basis set is also created from the information stored in the
    // checkpoint file
    // ...it needs to raid the checkpoint file to find those dimensions.
    // Create an SOBasis object using the basis set and integral factory.
    shared_ptr<SOBasisSet> soBasis(new SOBasisSet(aoBasis, integral));
    // Obtain block dimensions from the SO basis
    const Dimension dimension = soBasis->dimension();
    // The matrix factory can create matrices of the correct dimensions...
    shared_ptr<MatrixFactory> factory(new MatrixFactory);
    factory->init_with(dimension, dimension);
    // Form the one-electron integral objects from the integral factory
    shared_ptr<OneBodySOInt> sOBI(integral->so_overlap());
    shared_ptr<OneBodySOInt> tOBI(integral->so_kinetic());
    shared_ptr<OneBodySOInt> vOBI(integral->so_potential());
    // Form the one-electron integral matrices from the matrix factory
    SharedMatrix sMat(factory->create_matrix("Overlap"));
    SharedMatrix tMat(factory->create_matrix("Kinetic"));
    SharedMatrix vMat(factory->create_matrix("Potential"));
    SharedMatrix hMat(factory->create_matrix("One Electron Ints"));
    // Compute the one electron integrals, telling each object where to
    // store the result
    sOBI->compute(sMat);
    tOBI->compute(tMat);
    vOBI->compute(vMat);
    if(print > 5){
        sMat->print();
    }
    if(print > 3){
        tMat->print();
        vMat->print();
    }
    // Form h = T + V by first cloning T and then adding V
    hMat->copy(tMat);
    hMat->add(vMat);
    hMat->print();
    if(doTei){
        // 1. Obtain an object that knows how to compute two-electron AO
        // integrals.
        shared_ptr<TwoBodyAOInt> tb(integral->eri());
        // 2. Create an object that knows how to convert any two-body AO
        // integral to SO.
        shared_ptr<TwoBodySOInt> eri(new TwoBodySOInt(tb, integral));
        // 3. Find out how many SO shells we have.
        int nsoshell = soBasis->nshell();
        // 4. We to create an instance of our ERIPrinter
        ERIPrinter printer;
        // 5. Create an SOShellCombintationsIterator to step through the
        // necessary combinations
        SOShellCombinationsIterator shellIter(soBasis, soBasis, soBasis, soBasis);
        for (shellIter.first(); shellIter.is_done() == false; shellIter.next()) {
            // 6. Call the TwoBodySOInt object to compute integrals giving
            // it the
            // instance to our functor.
            eri->compute_shell(shellIter, printer);
        }
    }
    return Success;
}
}} // End namespaces

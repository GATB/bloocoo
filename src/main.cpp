/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

#include <DSK.hpp>
#include <Bloocoo.hpp>
#include <gatb/system/impl/System.hpp>
#include <gatb/tools/misc/impl/Property.hpp>

/********************************************************************************/

using namespace gatb::core;
using namespace std;

/********************************************************************************/

int main (int argc, char* argv[])
{
    // We define a try/catch block in case some method fails
    try
    {
        /** We create an instance of DSK class. */
        DSK dsk;

        /** We execute dsk. */
        dsk.run (argc, argv);

        /** We create an instance of Bloocoo class. */
        Bloocoo bloocoo;

        /** We execute dsk. */
        bloocoo.run (argc, argv);
    }

    catch (tools::misc::impl::OptionFailure& e)
    {
        if (e.getParser().saw("-h"))    {   e.getParser().displayHelp   (stdout);   }
        else                            {   e.getParser().displayErrors (stdout);   }
        return EXIT_FAILURE;
    }

    catch (system::Exception& e)
    {
        cerr << "EXCEPTION: " << e.getMessage() << endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

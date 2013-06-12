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
#include <gatb/tools/misc/impl/Tool.hpp>

/********************************************************************************/

using namespace gatb::core::tools;
using namespace std;

/********************************************************************************/

int main (int argc, char* argv[])
{
    // We define a try/catch block in case some method fails
    try
    {
        misc::impl::ToolComposite tool ("root");

        tool.add (new DSK());
        tool.add (new Bloocoo());

        tool.run (argc, argv);
    }

    catch (misc::impl::OptionFailure& e)
    {
        if (e.getParser().saw("-help")) {   e.getParser().displayHelp   (stdout);   }
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

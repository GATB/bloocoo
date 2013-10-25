/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

/********************************************************************************/

#include <Bloocoo.hpp>
#include <DSK.hpp>

using namespace std;

/********************************************************************************/

int main (int argc, char* argv[])
{
    // We define a try/catch block in case some method fails
    try
    {
        misc::impl::ToolComposite tool;

        tool.add (new DSK           ());
        tool.add (new Bloocoo       ());

        tool.run (argc, argv);
    }

    catch (misc::impl::OptionFailure& e)
    {
        e.getParser().displayErrors (stdout);
        e.getParser().displayHelp   (stdout);
        return EXIT_FAILURE;
    }

    catch (gatb::core::system::Exception& e)
    {
        cerr << "EXCEPTION: " << e.getMessage() << endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

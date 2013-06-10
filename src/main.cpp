/*****************************************************************************
 *   GATB : Genome Assembly Tool Box                                         *
 *   Authors: [R.Chikhi, G.Rizk, E.Drezen]                                   *
 *   Based on Minia, Authors: [R.Chikhi, G.Rizk], CeCILL license             *
 *   Copyright (c) INRIA, CeCILL license, 2013                               *
 *****************************************************************************/

#include <gatb/system/impl/System.hpp>
#include <gatb/tools/misc/impl/OptionsParser.hpp>
#include <gatb/tools/misc/impl/Property.hpp>

using namespace gatb::core::system::impl;
using namespace gatb::core::tools::misc;
using namespace gatb::core::tools::misc::impl;
using namespace std;

/********************************************************************************/

int main (int argc, char* argv[])
{
    OptionsParser parser;

    parser.add (new OptionOneParam (PROP_KMER_SIZE,   "size of a kmer",                       true));
    parser.add (new OptionOneParam (PROP_DB,          "URI of the bank",                      true));
    parser.add (new OptionOneParam (PROP_NB_CORES,    "number of cores",                      false));
    parser.add (new OptionOneParam (PROP_MAX_MEMORY,  "max memory",                           false));
    parser.add (new OptionOneParam (PROP_NKS,         "abundance threshold for solid kmers",  false));
    parser.add (new OptionOneParam (PROP_PREFIX,      "prefix URI for temporary files",       false));
    parser.add (new OptionNoParam  (PROP_QUIET,       "don't display exec information",       false));
    parser.add (new OptionOneParam (PROP_STATS_XML,   "dump exec info into a XML file",       false));
    
    // We define a try/catch block in case some method fails (bad filename for instance)
    try
    {
        /** We parse the command line arguments. */
        parser.parse (argc, argv);

        /** We get the options as a Properties object. */
        IProperties& props = parser.getProperties();

        /** We read properties from the init file (if any). */
        props.add (1, new Properties (System::info().getHomeDirectory() + string ("/.dskrc")));

#if 0
        /** We create an instance of DSK class. */
        DSK dsk (&props);

        /** We execute dsk. */
        dsk.execute ();

        /** We may have to dump execution information to stdout. */
        if (props[PROP_QUIET] == 0)
        {
            RawDumpPropertiesVisitor visit;
            dsk.getStats().accept     (&visit);
        }

        /** We may have to dump execution information to stdout. */
        if (props[PROP_STATS_XML] != 0)
        {
            XmlDumpPropertiesVisitor visit (props[PROP_STATS_XML]->getValue());
            dsk.getStats().accept     (&visit);
        }
#endif
    }

    catch (OptionFailure& e)
    {
        if (parser.saw("-h"))    {   parser.displayHelp   (stdout);   }
        else                     {   parser.displayErrors (stdout);   }
        return EXIT_FAILURE;
    }

    catch (gatb::core::system::Exception& e)
    {
        cerr << "EXCEPTION: " << e.getMessage() << endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}


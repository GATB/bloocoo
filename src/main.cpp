/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  INRIA
 *   Authors: G.Rizk, G.Benoit
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/

/********************************************************************************/

#include <Bloocoo.hpp>
//#include <DSK.hpp>

using namespace std;

void displayVersion(std::ostream& os){
	
	os << "* * * * * * * * * * * * * * * * * * * * * *" << endl;
	os << "* Bloocoo version "<< BLOOCOO_VERSION_MAJOR << "."
	<< BLOOCOO_VERSION_MINOR << "."
	<< BLOOCOO_VERSION_PATCH
	<< "                   *" << endl; //<< " AGPL licence" <<endl;
	os << "* Using gatb-core version "<< System::info().getVersion() <<  "           *" << endl;
	os << "* * * * * * * * * * * * * * * * * * * * * *" << endl;
}


/********************************************************************************/

int main (int argc, char* argv[])
{
	
	
	if(argc > 1 && (   strcmp(argv[1],STR_VERSION)==0 || strcmp(argv[1],"-v")==0    )     ){
		displayVersion(cout);
		return EXIT_FAILURE;
	}
	
	
    // We define a try/catch block in case some method fails
    try
    {
        Bloocoo().run (argc, argv);
    }

    catch (OptionFailure& e)
    {
        return e.displayErrors (cout);
    }

    catch (Exception& e)
    {
        cerr << "EXCEPTION: " << e.getMessage() << endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

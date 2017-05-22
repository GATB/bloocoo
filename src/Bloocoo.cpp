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

#include <Bloocoo.hpp>
//#include <DSK.hpp>  

#include <libgen.h>


// We use the required packages
using namespace std;

#define DEBUG(a)  printf a

//#define DEBION
/********************************************************************************/
// We define some string constants.
/********************************************************************************/
//const char* Bloocoo::STR_NB_MIN_VALID       = "-nbmin-valid";
const char* Bloocoo::STR_ION                = "-ion";
//const char* Bloocoo::STR_NB_VALIDATED_KMERS = "-nkmer-checked";
const char* Bloocoo::STR_ERR_TAB            = "-err-tab";
const char* Bloocoo::STR_MANUAL_SOLIDITY            = "-manual-solidity";

const char* Bloocoo::STR_RECALL            = "-high-recall";
const char* Bloocoo::STR_PRECISION            = "-high-precision";
const char* Bloocoo::STR_SLOW            = "-slow";
const char* Bloocoo::STR_MAX_TRIM            = "-max-trim";

const char* Bloocoo::STR_HISTO_ONLY            = "-count-only";
const char* Bloocoo::STR_FROM_H5            = "-from-h5";
const char* Bloocoo::STR_BIT_BLOOM_PER_KMER = "-nbits-bloom";


// these 2 tabs are now known globally
//char bin2NT[4] = {'A','C','T','G'};
//char binrev[4]    = {2,3,0,1};

//char bin2NTrev[4] = {'T','G','A','C'};



/********************************************************************************/
//fonctions for correction
/********************************************************************************/



//#define PRINT_DEBUG
//#define PRINT_STATS // only in single thread mode

#define PRINT_LOG_MULTI 0


//#define SERIAL // to force serial mode





int NT2int(char nt)
{
    int i;
    i = nt;
    i = (i>>1)&3; // that's quite clever, guillaume.
    return i;
}


inline unsigned int popcount_32(unsigned int x)
{
    unsigned int m1 = 0x55555555;
    unsigned int m2 = 0x33333333;
    unsigned int m4 = 0x0f0f0f0f;
    unsigned int h01 = 0x01010101;
    x -= (x >> 1) & m1;               /* put count of each 2 bits into those 2 bits */
    x = (x & m2) + ((x >> 2) & m2);   /* put count of each 4 bits in */
    x = (x + (x >> 4)) & m4;          /* put count of each 8 bits in partie droite  4bit piece*/
    return (x * h01) >> 24;           /* returns left 8 bits of x + (x<<8) + ... */
}


inline unsigned int popcount_64(uint64_t x)
{
    
    unsigned int low = x & 0xffffffff ;
    unsigned int high = ( x >> 32) & 0xffffffff ;
    
    return (popcount_32(low) + popcount_32(high));
}



//return 1 if homopo at first pos,
//2 if homopo at second pos
//detect homopo of at least 2nt (todo configure this)
int contains_homopolymer (kmer_type kmer, int sizeKmer)
{
    
    kmer_type one; one.setVal(1);
    kmer_type  kmerMask=(one<<(sizeKmer*2))-1;
    
    // kmer_type rfin = 0xffffffffffffffc0; //elim 3 derniers
    kmer_type rfin; rfin.setVal( 0xffffffffffffff00); //elim 4 derniers
    
    kmer_type m; m.setVal(0x5555555555555555);
    
    //    char kmerbuff[sizeKmer+1];
    //    code2seq(kmer,kmerbuff);
    //    printf("kmer: %s \n",kmerbuff);
    
    //  std::cout << kmer.toString(sizeKmer)  << std::endl; ;
    
    
    kmer_type x = kmer;
    //printf("%llx\n",x.getVal());
    
    x = ~ (kmer ^ (kmer << 2)) ;
    //printf("%llx   (nxor d2 )\n",x.getVal());
    
//    //pour 3 nt au moins
//    kmer_type y = ~ (kmer ^ (kmer << 4)) ; //dec de 2 nt
//    x= x & y ;
//    /////
    
//    //pour 4 nt au moins
//    kmer_type z = ~ (kmer ^ (kmer << 6)) ; //dec de 3 nt
//    x= x & z ;
    
    ///////////
    //printf("%llx   (nxor d246 )\n",x.getVal());
    
    
    x = (x & (x >> 1)) & m ;
    //printf("%llx   (fin )\n",x.getVal());
    
    kmer_type three; three.setVal(3);
    kmer_type  maskFirstNt=(three<<((sizeKmer-1)*2)); // select first nt
    kmer_type  maskSecondNt=(three<<((sizeKmer-2)*2));
    
    if((x & maskFirstNt).getVal()  >0) return 1;
    if((x & maskSecondNt).getVal() >0 ) return 2;
    
    return 0 ;
    //x= (x & kmerMask ) & rfin;
    //puis popcount
    //printf("%llx   (mm )\n--\n",x.getVal());
    
    
    // printf(" popcnt %i\n",popcount_64(x ));
    
    return (popcount_64(x.getVal() )>0);
}

/*
//default constructor
CorrectReads ()
: _bloom(NULL), _outbank(NULL), _bloocoo(NULL),
_total_nb_errors_corrected (NULL), _local_nb_errors_corrected(0),
_total_nb_ins_corrected (NULL), _local_nb_ins_corrected(0),
_total_nb_del_corrected (NULL), _local_nb_del_corrected(0),
_bankwriter(NULL),_temp_nb_seq_done(0)
{

  //   printf("------- CorrectReads default Constructor  %p --------- tid %i  \n",this,_thread_id);
	_tab_multivote = (unsigned char *) malloc(TAB_MULTIVOTE_SIZE*sizeof(unsigned char)); // pair of muta  = 16 nt *128 pos * 16 (max dist)
	_newSeq = 0;
	if(_bloocoo->_ion_mode)
	{
	_newSeq = new Sequence ();
	}
}*/











/********************************************************************************/
/* Class Bloocoo for read correction : takes as input a bag of solid kmers (dsk result),
 insert that in a bloom filter, and correct read form it*/
/********************************************************************************/
Bloocoo::Bloocoo () : Tool("bloocoo"), _kmerSize(27), _inputBank (0), errtab_mutex(0)
{
    // _seq_num = 0;
    
    #ifdef PRINT_STATS
		for(int i=0; i<NB_CORRECTION_METHODS; i++){
			__correction_methods_successes[i] = 0;
            __correction_methods_calls[i] = 0;
            
		}
	#endif

    /** We add the sorting count options and hide all of them by default and display one some of them. */
    getParser()->push_back (SortingCountAlgorithm<>::getOptionsParser(), 1);

    if (IOptionsParser* input = getParser()->getParser (STR_URI_INPUT))  {  input->setName (STR_URI_FILE);  }

    /** We add options specific to this tool. */
    getParser()->push_back(new OptionNoParam (Bloocoo::STR_RECALL, "correct a lot but can introduce more mistakes", false));
    getParser()->push_back(new OptionNoParam (Bloocoo::STR_PRECISION, "correct safely, correct less but introduce less mistakes", false));
    getParser()->push_back(new OptionNoParam (Bloocoo::STR_SLOW, "slower but better correction", false));
    //getParser()->push_back (new OptionOneParam (Bloocoo::STR_NB_MIN_VALID,        "min number of kmers to valid a correction", false,"3"));
    //getParser()->push_back (new OptionOneParam (Bloocoo::STR_NB_VALIDATED_KMERS,  "number of kmers checked when correcting an error", false,"4"));
    getParser()->push_back (new OptionNoParam (Bloocoo::STR_ION, "ion correction mode", false));
    getParser()->push_back (new OptionNoParam (Bloocoo::STR_ERR_TAB, "generate error tabs", false));
    getParser()->push_back(new OptionOneParam(STR_MAX_TRIM, "max number of bases that can be trimmed per read", false));

	getParser()->push_back(new OptionOneParam(STR_BIT_BLOOM_PER_KMER, "(advanced option) nb bits per kmer for bloom filter", false,"12")); //

	
	
	
	//had to add this becasue I cannot see if STR_KMER_ABUNDANCE_MIN was set or not (getParser()->saw() always return true because it has a default value)
	//and I cannot change STR_KMER_ABUNDANCE_MIN default value because it is in gatb-core
	getParser()->push_back (new OptionNoParam (Bloocoo::STR_MANUAL_SOLIDITY, "switch to manual solidity (use -abundance-min to set it)", false));

	
	getParser()->push_back(new OptionNoParam (Bloocoo::STR_FROM_H5, "do not re-compute kmer counts, suppose h5 file already computed (with previous run with -count-only)", false));
	getParser()->push_back(new OptionNoParam (Bloocoo::STR_HISTO_ONLY, "do not correct, only count kmers", false));

	
}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
Bloocoo::~Bloocoo ()
{
    /** Cleanup. */
    if (errtab_mutex != 0)  { delete errtab_mutex; }
}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
void Bloocoo::execute ()
{
	cout << "Executing Bloocoo" << endl;
	//__error_detected = 0;
	
    /*************************************************/
    // We set some attributes (shortcuts).
    /*************************************************/
    _kmerSize           = getInput()->getInt (STR_KMER_SIZE);
    chooseCorrectionParams();
    //_nb_kmers_checked   = getInput()->getInt (STR_NB_VALIDATED_KMERS);
    //_nb_min_valid = getInput()->getInt (STR_NB_MIN_VALID);

    string inputFilename = getInput()->getStr(STR_URI_FILE);

    _solidFile = getInput()->get(STR_URI_OUTPUT) ?
        getInput()->getStr(STR_URI_OUTPUT) + ".h5"  :
        System::file().getBaseName (inputFilename) + ".h5"; //_inputFilename instead of prefix GR

	if( _fromH5Mode )
	{
		FILE* testfile = fopen (_solidFile.c_str(), "r");
		if (testfile == 0)
		{
			fprintf(stderr,"File %s not found, cannot start from h5 file\n",_solidFile.c_str());
			exit(1);
		}
		else
		{
			fclose( testfile);
		}
	}
	
    /*************************************************/
    // Sorting count part
    /*************************************************/
	if(! _fromH5Mode) //in from h5 mode, suppose the h5 file is present
    {
        /** We open the input bank. */
        IBank* inputBank = Bank::open(inputFilename);
        LOCAL (inputBank);

        /** We create a SortingCountAlgorithm instance. */
        SortingCountAlgorithm<> sortingCount (inputBank, getInput());

        sortingCount.getInput()->add (0, STR_VERBOSE, getInput()->getStr(STR_VERBOSE));

        /** We execute the sorting count. */
        sortingCount.execute();
		
		
		if(_countOnlyMode)
		{
			printf("Count only mode : kmer counts and histogram of occurences are in the  %s file \n",_solidFile.c_str());
			exit(0);
		}
    }
    
    /*************************************************/
    /** We create a bloom with inserted solid kmers. */
    /*************************************************/
    _bloom = createBloom ();
    LOCAL (_bloom);
    
    //iterate over initial file
    	IBank* inbank = Bank::open(inputFilename);

    //cout <<  inputFilename << std::endl;
    
    _max_trim = 0;
    if(getInput()->get(STR_MAX_TRIM)){
		_max_trim = getInput()->getInt(STR_MAX_TRIM);
	}
	
	//cout << "Max trim: " << _max_trim << endl;
    
    _ion_mode= false;
    if ( getParser()->saw (Bloocoo::STR_ION) )
    {
        _ion_mode = true;
    }
    
    
    
    bool fastq_mode= false;
    if( inputFilename.find("fastq") !=  string::npos )  fastq_mode =true;
    //printf("-- %s -- fq mode %i \n",getInput()->getStr(DSK::STR_URI_DATABASE).c_str(),fastq_mode);
        if( inputFilename.find("fq") !=  string::npos )  fastq_mode =true;

	printf("fastq mode %i \n",fastq_mode);
	
    /*************************************************/
    // We create a sequence iterator for the bank
    /*************************************************/
	Iterator<Sequence>* itSeq =
	//pourquoi ce create iterator fait perdre la compo de literator ?  ensuite  itSeq->getComposition() renvoit 1 truc seul : il y avait bug gatb-core
	createIterator<Sequence> (
                                                          inbank->iterator(),
                                                          inbank->estimateNbItems(),
                                                          "Iterating and correcting sequences"
                                                          );
	
    LOCAL (itSeq);
    
    /*************************************************/
    // We create the corrected file
    /*************************************************/
	string prefix;
    string outputFilename;
    if(getInput()->get(STR_URI_OUTPUT)){
		outputFilename = getInput()->getStr(STR_URI_OUTPUT);
		prefix = System::file().getBaseName(outputFilename); // FIXME or delete this comment: won't this override the path? e.g. if the user set an output file to another folder 
	}
	else{
		/** We get the basename from the provided URI (ie remove directory path and suffix). */
		//<<<<<<< Updated upstream
		prefix = System::file().getBaseName(inputFilename);
		//=======
		//    string prefix = System::file().getBaseName (getInput()->getStr(DSK::STR_URI_DATABASE));
		//
		//>>>>>>> Stashed changes
		/** We set the filename as the base name + a specific suffix. */
		//string fileName;
		if(fastq_mode)
			outputFilename = prefix + string("_corrected.fastq");
		else
			outputFilename = prefix + string("_corrected.fasta");
	}


    
    bool gz_mode= false;
	
    
    //BagCache<Sequence> *  outbank_cache = new BagCache<Sequence> (&outbank,10000);
    
    u_int64_t total_nb_errors_corrected = 0;
    u_int64_t total_nb_ins_corrected = 0;
    u_int64_t total_nb_del_corrected = 0;

    int nb_corrector_threads_living;
    
    _wantErrTabFile = getInput()->get(Bloocoo::STR_ERR_TAB) != 0;
    if(_wantErrTabFile){
		
        //file with list of errors, for testing purposes
        string ferrfile = prefix + string ("_bloocoo_corr_errs.tab");

        _errfile = System::file().newFile (ferrfile, "wb");

        //file with list of errors, for testing purposes, with name of method responsible of correction
        string ferrfile2 = prefix + string ("_bloocoo_corr_errs_full.tab");

        _errfile_full = System::file().newFile (ferrfile2, "wb");
    }
    
	#if PRINT_LOG_MULTI
		string ferrfile3 = prefix + string ("_debug.tab");
		_debug = System::file().newFile (ferrfile3, "wb");
	#endif
	
	
	std::vector<string> outbanknames;
    /*************************************************/
    // We iterate over sequences and correct them
    /*************************************************/
    {
        TIME_INFO (getTimeInfo(), "sequences correction");
        
        //   CorrectReads fct2(*bloom, outbank, *this,total_nb_errors_corrected) ;
        //   itSeq->iterate (fct2); //
        //
        
        nb_corrector_threads_living = 0 ;
        
		#ifdef SERIAL
			setDispatcher (new SerialDispatcher());
		#else
			setDispatcher (  new Dispatcher (getInput()->getInt(STR_NB_CORES)) );
		#endif
		
		
		std::vector<Iterator<Sequence>*> itBanks =  itSeq->getComposition();

		//printf("itseq nb compo %i \n",itBanks.size());
		
		for (size_t i=0; i<itBanks.size(); i++)
		{
		//	printf("should iterate bank i %i   bankalbum id  %s ,sub bank id %s  \n",i, inbank->getId().c_str(),inbank->getIdNb(i).c_str() );
			
			

			fastq_mode =false;
			string bankname;
			if(itBanks.size()>1)
				bankname = (inbank->getIdNb(i)); // System::file().getBaseName
			else
				bankname = inbank->getId(); // System::file().getBaseName(inbank->getId()); // was System::file().getBaseName
			
			if( bankname.find("fastq") !=  string::npos )  fastq_mode =true;
			if( bankname.find("fq") !=  string::npos )     fastq_mode =true;
			
			
			printf("bank name %s fastq mode %i \n",bankname.c_str(),fastq_mode);
			
			////name management
			if(!(getInput()->get(STR_URI_OUTPUT)))
			{
				if(itBanks.size()>1)
					prefix = System::file().getBaseName(inbank->getIdNb(i));
				else
					prefix = System::file().getBaseName(inbank->getId());
			}
			else //name provided by user
			{
				if(itBanks.size()>1)
					prefix =  System::file().getBaseName(getInput()->getStr(STR_URI_OUTPUT)) + Stringify::format ("_%i_", i);
				else
					prefix =  System::file().getBaseName(getInput()->getStr(STR_URI_OUTPUT));
			}
			
			if(  !(getInput()->get(STR_URI_OUTPUT)) )
			{
				if(fastq_mode)
					outputFilename = prefix + string("_corrected.fastq");
				else
					outputFilename = prefix + string("_corrected.fasta");
			}
			else //name provided by user
			{
				if(itBanks.size()>1)
				{
					if(fastq_mode)
						outputFilename = prefix + string(".fastq");
					else
						outputFilename = prefix + string(".fasta");
				}
				else
					outputFilename =  getInput()->getStr(STR_URI_OUTPUT);
			}
			/////////
			
			outbanknames.push_back(outputFilename);
			
			if(itBanks.size()>1)
				cout << "\tInput filename: " << inbank->getIdNb(i) << endl;
			else
				cout << "\tInput filename: " << inbank->getId() << endl;
			
			cout << "\tOutput filename: " << outputFilename << endl;
			
			
			
			char msg[1000];
			if(itBanks.size()>1)
				snprintf(msg,1000,"Correcting sequences, file %i",i+1);
			else
				snprintf(msg,1000,"Iterating and correcting sequences");

			
			
		 	itSeq =
					createIterator<Sequence> (
									itBanks[i],
									inbank->estimateNbItemsBanki(i),
									msg
									  );
			
			
			
			
			BankFasta outbank (outputFilename, fastq_mode,gz_mode);

			
			getDispatcher()->iterate (itSeq,  CorrectReads (_bloom, &outbank , this, &total_nb_errors_corrected, &total_nb_ins_corrected, &total_nb_del_corrected,getInput()->getInt(STR_NB_CORES),&nb_corrector_threads_living),10000); // was 10000
			
			itBanks[i]->finalize();
			
		}
		
        //replace outbank_cache with &outbank to remove usage of bagcache
 // , 10000
        
        // printf("---after dispatcher iterate ---\n");
    }

    /*************************************************/
    // We gather some statistics.
    /*************************************************/
    getInfo()->add (1, "result");
    getInfo()->add (2, "nb errors corrected", "%ld", total_nb_errors_corrected);
    if(_ion_mode)
    {
        getInfo()->add (2, "nb ins corrected", "%ld", total_nb_ins_corrected);
        getInfo()->add (2, "nb del corrected", "%ld", total_nb_del_corrected);
    }
	string outfile ;

	for (size_t i=0; i<inbank->getCompositionNb(); i++)
	{
		if(i>0) outfile+= " , ";
		outfile+=  outbanknames[i];
	}
	
	if(inbank->getCompositionNb()>1)
    	getInfo()->add (2, "corrected file", outfile);
	else
		getInfo()->add (2, "corrected file", outputFilename);
	
    #ifdef PRINT_STATS
    	printf("Correction methods stats (TwoSided, Agressive, Vote, MultiMutateVote):\n");
		for(int i=0; i<NB_CORRECTION_METHODS; i++){
			printf("%i    /  %i \n", __correction_methods_successes[i], __correction_methods_calls[i] );
		}
	#endif

}

/*********************************************************************
 ** secure mode: correction with high precision
 ** slow mode: more correction loop for each read
 ** _nb_less_restrictive_correction: perform n correction loop while the read is not fully corrected (n is the value of _nb_less_restrictive_correction)
 ** _only_decrease_nb_min_valid: If true, each correction loop is performed with a reduced _nb_min_valid (meaning a better recall).
 **                              If false, both _nb_min_valid and _nb_kmers_checked are reduced (meaning a better precision)
 *********************************************************************/
void Bloocoo::chooseCorrectionParams(){
    bool isRecallMode = getParser()->saw(Bloocoo::STR_RECALL);
    bool isPrecisionMode = getParser()->saw(Bloocoo::STR_PRECISION);
    bool isSlowMode = getParser()->saw(Bloocoo::STR_SLOW);
	
	 _fromH5Mode = getParser()->saw(Bloocoo::STR_FROM_H5);
	 _countOnlyMode = getParser()->saw(Bloocoo::STR_HISTO_ONLY);
	_manualSolidity = getParser()->saw(Bloocoo::STR_MANUAL_SOLIDITY);
	
	
	if(isSlowMode){
		_nb_less_restrictive_correction = 2;
		if(isRecallMode){ //8-4 8-2
			_nb_kmers_checked = 8;
			_nb_min_valid = 4;
			_only_decrease_nb_min_valid = true;
		}
		else if(isPrecisionMode){ //8-7 6-5
			_nb_kmers_checked = 8;
			_nb_min_valid = 7;
			_only_decrease_nb_min_valid = false;
		} 
		else{ //8-7 8-5
			_nb_kmers_checked = 8;
			_nb_min_valid = 6;
			_only_decrease_nb_min_valid = true;
		}
	}
	else{
		_nb_less_restrictive_correction = 1;
		_only_decrease_nb_min_valid = false;
		if(isRecallMode){ 
			_nb_kmers_checked = 5;
			_nb_min_valid = 2;
		}
		else if(isPrecisionMode){ 
			_nb_kmers_checked = 5;
			_nb_min_valid = 5;
		} 
		else{
			_nb_kmers_checked = 5;
			_nb_min_valid = 4;
		}
	}
    

	_max_multimutation_distance = 6;
	
	cout << "Correction params are:" << endl;
	if(isRecallMode)
		cout << "\tMode: High-recall" << endl;
	else if(isPrecisionMode)
		cout << "\tMode: High-precision" << endl;
	else
		cout << "\tMode: Default" << endl;
	if(isSlowMode)
		cout << "\tSlow mode: Yes" << endl;
	else
		cout << "\tSlow mode: No" << endl;
	cout << "\tKmers checked: " << _nb_kmers_checked << endl;
	cout << "\tMin kmers valid: " << _nb_min_valid << endl;
	cout << "\tDecrease only kmers min valid: " << _only_decrease_nb_min_valid << endl;
	cout << "\tCorrection loop: " << _nb_less_restrictive_correction << endl;
}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
IBloom<kmer_type>* Bloocoo::createBloom ()
{
    TIME_INFO (getTimeInfo(), "fill bloom filter");
    
    double lg2 = log(2);
    float NBITS_PER_KMER = log (16*_kmerSize*(lg2*lg2))/(lg2*lg2);
    NBITS_PER_KMER = 12;

	if(getInput()->get(STR_BIT_BLOOM_PER_KMER)){
		NBITS_PER_KMER = getInput()->getInt(STR_BIT_BLOOM_PER_KMER);
	}
	printf("NBITS per kmer %f \n",NBITS_PER_KMER);

	
	int auto_cutoff = 0 ;
	u_int64_t nbs = 0 ;
	u_int64_t nb_kmers_infile;
	int _nks;
	
    /** We retrieve the solid kmers from the storage content (likely built by DSK). */
    Storage* storage = StorageFactory(STORAGE_HDF5).create (_solidFile, false, false);
    LOCAL (storage);

    /** We get an iterator over all the [kmer,abundance]  (ie. from all solid partitions). */
    Partition<kmer_count> & solidCollection = storage->root().getGroup("dsk").getPartition<kmer_count> ("solid");

    /** We get the number of solid kmers. */
    u_int64_t solidFileSize = solidCollection.getNbItems();
	nb_kmers_infile = solidCollection.getNbItems();

	

		//	if( ! getParser()->saw(STR_KMER_ABUNDANCE_MIN)){ // this does not work ! always return true because it has a default value
		if(!_manualSolidity){
		
		//retrieve cutoff
		//printf("retrieving cutoff \n");
		
		Collection<NativeInt64>& cutoff  = storage->getGroup("histogram").getCollection<NativeInt64> ("cutoff");
		Iterator<NativeInt64>* iter = cutoff.iterator();
		LOCAL (iter);
		for (iter->first(); !iter->isDone(); iter->next())  {
			auto_cutoff = iter->item().toInt();
		}
		
		//retrieve nb solids
		
		Collection<NativeInt64>& storagesolid  = storage->getGroup("histogram").getCollection<NativeInt64> ("nbsolidsforcutoff");
		Iterator<NativeInt64>* iter2 = storagesolid.iterator();
		LOCAL (iter2);
		for (iter2->first(); !iter2->isDone(); iter2->next())  {
			nbs = iter2->item().toInt();
		}
		
		_nks = auto_cutoff;
	}
	else
	{
		auto_cutoff =0;
		_nks =  getInput()->getInt(STR_KMER_ABUNDANCE_MIN);
		nbs  = nb_kmers_infile;
	}
	
	
	
    u_int64_t estimatedBloomSize = (u_int64_t) ((double)nbs * NBITS_PER_KMER);
    if (estimatedBloomSize ==0 ) { estimatedBloomSize = 1000; }
    
    /** We create the kmers iterator from the solid file. */
    Iterator<kmer_count>* itKmers = createIterator<kmer_count> (
        solidCollection.iterator(),
        solidFileSize,
        "fill bloom filter"
    );
    LOCAL (itKmers);
	
	
	if(auto_cutoff)
		cout << "Abundance threshold: " << auto_cutoff << " (auto)    (nb solid kmers: " << nbs << ")"<< endl;
	else
		cout << "Abundance threshold: " << _nks << " (user set)  (nb solid kmers: " << nbs << ")"<< endl;
	
	
	
    /** We instantiate the bloom object. */
    IProperties* stats = new Properties();  LOCAL (stats);
    BloomBuilder<> builder (estimatedBloomSize, 7, _kmerSize, tools::misc::BLOOM_CACHE,getInput()->getInt(STR_NB_CORES),auto_cutoff);
    IBloom<kmer_type>* bloom = builder.build (itKmers, stats);
    
    getInfo()->add (1, "bloom");
    getInfo()->add (2, stats);

    /** We return the created bloom filter. */
    return bloom;
}


/********************************************************************************/
/* Class CorrectReads used to correct a single read. Bloocoo creates CorrectReads instances
 * to perform the correction in multiple threads. Each instances of this classes
 * is synchronized*/
/********************************************************************************/
CorrectReads::CorrectReads(IBloom<kmer_type>* bloom, Bag<Sequence>* outbank, Bloocoo * bloocoo, u_int64_t* nb_errors_corrected, u_int64_t* nb_ins_corrected, u_int64_t* nb_del_corrected, int nb_cores, int * nbliving)
: _bloom(bloom), _outbank(outbank), _bloocoo(bloocoo),
_total_nb_errors_corrected (nb_errors_corrected),_total_nb_ins_corrected (nb_ins_corrected),
_total_nb_del_corrected (nb_del_corrected),_local_nb_errors_corrected(0),
_synchro(0),_temp_nb_seq_done(0), _nb_living(nbliving), _bankwriter(0)
{
	setBankWriter (new OrderedBankWriter(outbank,nb_cores*10000));// was 10000
	_thread_id = __sync_fetch_and_add (_nb_living, 1);
	_local_nb_ins_corrected =0;
	_local_nb_del_corrected =0;
	_tab_multivote = (unsigned char *) malloc(TAB_MULTIVOTE_SIZE*sizeof(unsigned char)); // pair of muta  = 16 nt *128 pos * 16 (max dist)
     // printf("------- CorrectReads Custom Constructor  %p --------- tid %i \n",this,_thread_id);
	
	setSynchro (System::thread().newSynchronizer());

	_nb_less_restrictive_correction = _bloocoo->_nb_less_restrictive_correction;
    _max_multimutation_distance = _bloocoo->_max_multimutation_distance;
    _only_decrease_nb_min_valid = _bloocoo->_only_decrease_nb_min_valid; // false = descend les 2, recall plus faible
	_kmerSize = _bloocoo->_kmerSize;
	_nb_kmers_checked = _bloocoo->_nb_kmers_checked;
	_nb_min_valid = _bloocoo->_nb_min_valid;
    _wantErrTabFile = _bloocoo->_wantErrTabFile;
    
	_newSeq =0;
	if(_bloocoo->_ion_mode)
	{
		//_tempseq = (char *) malloc(sizeof(char)*10000); // todo chande max seq size ..
		_newSeq = new Sequence (); //readseq

	}
	
	//_hashAnchorKmers = new OAHash<kmer_type>(1000000);
	
}

//copy construct
CorrectReads::CorrectReads(const CorrectReads& cr) //called by dispatcher iterate to create N functors
    : _synchro(0), _bankwriter(0)
{
	//functors share smee bloom, bankwriter, bloocoo and synchronizer
	_bloom = cr._bloom;
	_outbank = cr._outbank;
	_bloocoo = cr._bloocoo;
	_nb_living = cr._nb_living;

	/** We share the synchronizer and the bank writer (both SmartPointers) between the different instances of CorrectReads. */
	setSynchro    (cr._synchro);
    setBankWriter (cr._bankwriter);

	_newSeq = 0;
	_tab_multivote = (unsigned char *) malloc(TAB_MULTIVOTE_SIZE*sizeof(unsigned char)); // pair of muta  = 16 nt *128 pos * 16 (max dist)
	
	_total_nb_errors_corrected = cr._total_nb_errors_corrected;
	_total_nb_ins_corrected = cr._total_nb_ins_corrected;
	_total_nb_del_corrected = cr._total_nb_del_corrected;

	_local_nb_errors_corrected =0;
	_local_nb_ins_corrected =0;
	_local_nb_del_corrected =0;

	_temp_nb_seq_done = 0;
	_thread_id = __sync_fetch_and_add (_nb_living, 1);
	// printf("------- CorrectReads copy construct  from  %p  into %p --------- tid %i \n",&cr,this,_thread_id);
	
	_nb_less_restrictive_correction = cr._nb_less_restrictive_correction;
	_max_multimutation_distance = cr._max_multimutation_distance;
	_only_decrease_nb_min_valid = cr._only_decrease_nb_min_valid;
	_kmerSize = cr._kmerSize;
	_nb_kmers_checked = cr._nb_kmers_checked;
	_nb_min_valid = cr._nb_min_valid;
	_wantErrTabFile = cr._wantErrTabFile;
	
	if(_bloocoo->_ion_mode)
	{
		_newSeq = new Sequence ();
	}
}
//assign
//    CorrectReads& operator=(const CorrectReads& cr)
//    {
//     //   printf("------- CorrectReads assign operator from  %p   ---------\n",&cr);
//        return *this;
//    }

CorrectReads::~CorrectReads ()
{
 getSynchro()->lock()  ;
	/** We increase the global number of corrected errors. */
	//  printf("local nb errors %lli \n",_local_nb_errors_corrected);
	 __sync_fetch_and_add (_total_nb_errors_corrected, _local_nb_errors_corrected);
	 __sync_fetch_and_add (_total_nb_ins_corrected, _local_nb_ins_corrected);
	 __sync_fetch_and_add (_total_nb_del_corrected, _local_nb_del_corrected);
	
   // printf("tot nb ins c %lli loc %lli \n",*_total_nb_ins_corrected,_local_nb_ins_corrected);

	
	#ifndef SERIAL
		_bankwriter->incDone(_temp_nb_seq_done);
		_bankwriter->waitForWriter();
	#endif
	
	
	int nb_remaining = __sync_fetch_and_add (_nb_living, -1);

	//printf("------- CorrectReads Destructor  %p --------- nb_remaining  %i \n",this,nb_remaining);

	
	
	#ifndef SERIAL
		if(nb_remaining==1)
		{
		//	printf("should flush writer  nbr %i \n",nb_remaining);
			_bankwriter->FlushWriter();
		}
	#endif
	
	if (_tab_multivote) { free(_tab_multivote); } // pourquoi plante ? lobjet correct read est il copié qq part ?
	if (_newSeq)        { delete _newSeq;       }
 getSynchro()->unlock()  ;

	/** We release one token of the smart pointers. */
	setSynchro    (0);
    setBankWriter (0); //peut etre ceci est fait avant que le writer thread ne soit killé ?
}

//no const otherwise error with tKmer.setData
void CorrectReads::operator()(Sequence& sequence){
	_sequence = &sequence;
	
	_readseq = _sequence->getDataBuffer(); // the nucleotide sequence of the read
	_readSize = _sequence->getDataSize();
	_kmerCount = _readSize - _kmerSize + 1;
	
	if(_readSize <= _kmerSize){
		writeSequence();
		return;
	}
	
	_corrected_pos.clear();
	_kmers.clear();
	_nb_kmer_offset = -2;
	_use_newSeq = false;

	if(_bloocoo->_ion_mode){
		executeIonCorrection();
	}
	else{
		execute2();
	}

	writeSequence();
}

/*
//Not used anymore
void CorrectReads::execute(){
	
	int nb_tries =0;
	
	
	
	int nb_checked;
	KmerModel model (_kmerSize,KMER_DIRECT);
	KmerModel::Iterator itKmer (model);
	kmer_type current_kmer;
	kmer_type current_kmer_min;
	
	_pos_homopo = 0;
	_last_wrong_kmer = 0;
	
	for(int j=0; j<_nb_less_restrictive_correction; j++){
		_nb_kmer_offset += 2;
		
		
		int nb_kmer_trusted = 0;
		_continue_correction = true;
		bool first_gap = true;
		int ii=0;
		
		
		while(_continue_correction && nb_tries <5)
		{
			nb_tries ++;
			nb_kmer_trusted = 0;
			_continue_correction = false;
			first_gap = true;
			ii = 0;
			
			if(_bloocoo->_ion_mode){
				itKmer.setData (_sequence->getData());
				_readSize = _sequence->getDataSize();
				_pos_homopo = 0;
			   // printf("current seq size %i \n",_readSize);
			}
			else{
				itKmer.setData (_sequence->getData());
			}
			#ifdef DEBION
				printf("%.*s   %i  \n",(int)_sequence->getDataSize(),_sequence->getDataBuffer(),_sequence->getIndex());
			#endif
			int untrusted_zone_size = 0;
			int trusted_zone_size = 0;
			int real_untrusted_zone_size = 0;
			
			_kmers.clear();
			
			
			//if(PRINT_DEBUG){ _bloocoo->print_read_correction_state(&model, _sequence);}
			
			// getSynchro()->lock()  ;
			
			// We iterate the kmers of this sequence
			for (itKmer.first(); !itKmer.isDone(); itKmer.next(),ii++)
			{
				
				int nb_errors_cor = 0;
				
				current_kmer = *itKmer;
				//kmers[ii] = &(*itKmer);
				_kmers.push_back(current_kmer);
				
				current_kmer_min = min(revcomp(current_kmer, _kmerSize), current_kmer);
				
				
				
				//Si on veut checker le dernier kmer du read et qu'on est dans une zone untrusted alors on considere
				//automatiquement ce dernier kmer comme untrusted et on laisse les démarches de correction de fin de read s'effectuer.
				bool is_last_kmer_indexed_after_hole = (ii==_readSize-_kmerSize && trusted_zone_size==0);
				
				
				if (!is_last_kmer_indexed_after_hole && _bloom->contains(current_kmer_min)) //kmer is solid
				{
					
					nb_kmer_trusted += 1;
					
					trusted_zone_size += 1;
					
					//beginning of indexed zone
					if(trusted_zone_size == 2)
					{
						if(_bloocoo->_ion_mode)
						{
							ionCorrection(ii);
							if(_continue_correction) break;
						}
						if (untrusted_zone_size>1){
							
							//_bloocoo->__error_detected += 1;
							//if(ii-2 == 30){
							//	printf("%i %i %i %i %i\n", (untrusted_zone_size-1) == _kmerSize, ii > _kmerSize, ii-2, untrusted_zone_size, _kmerSize);
							//}
							
							//two-sided conservative correction
							// this should be an isolated error, middle of the read
							// if error at pos _kmerSizecorrected_pos+1+ii, could work in theory but would need kmer_begin+1
							if(!_bloocoo->_ion_mode)
							{
								if(((untrusted_zone_size-1) == _kmerSize)  &&  (ii > _kmerSize+2)){
									
									nb_errors_cor = twoSidedCorrection(ii-2);
									
								}
								
								if(nb_errors_cor == 0){
									nb_checked = _nb_kmers_checked;
									nb_checked = min(_nb_kmers_checked, untrusted_zone_size-1); // -2 changed to -1 : was caus of pb for sides
									nb_checked = max(nb_checked, 0);
									
									//if first gap, cannot correct from the left side of the gap
									//2eme condition: Agressive right peut corriger le debut si le trou est plus long que _kmerSize mais est-ce utile?
									if(!first_gap)// || (first_gap && untrusted_zone_size > _kmerSize))
									{
										nb_errors_cor = aggressiveCorrection(ii-untrusted_zone_size, nb_checked, RIGHT);
									}
									
									
									if(nb_errors_cor == 0){
										nb_errors_cor = aggressiveCorrection(ii-2, nb_checked, LEFT);
									}
									
									if(nb_errors_cor == 0){
										nb_errors_cor = voteCorrectionInUntrustedZone(ii- untrusted_zone_size, ii-2, nb_checked);
									}
									
									if(nb_errors_cor == 0){
										nb_errors_cor = multiMutateVoteCorrectionInUntrustedZone(ii- untrusted_zone_size, ii-2, nb_checked,-1,ii-2);
									}
									
									
								}
							}
						}
						
					}
					
					if(trusted_zone_size==2){
						first_gap = false;
					}
					
					////////////traitement d'un faux positif du bloom isolé
					//   if(trusted_zone_size==1) tai_previous_break = untrusted_zone_size; // begin of indexed zone    ( if untrusted_zone_size ?)
					if (trusted_zone_size > 1) // do not reset  untrusted_zone_size if a single positive (might be a FP)
					{
						untrusted_zone_size = 0; //reset previous hole size
					}
					if (trusted_zone_size==1)// begin of indexed zone
					{
						untrusted_zone_size++; //  treat a single solidkmer  as an erroneous kmer
					}
					//////////
					
					//version sans traitement du FP du bloom
					//if(trusted_zone_size==1){
					//	untrusted_zone_size = 0;
					//}
					
				}
				else //kmer contains an error
				{
					untrusted_zone_size += 1;
					trusted_zone_size = 0;
					
					
					//todo put this in outside func
					if(_bloocoo->_ion_mode)
					{
						_last_wrong_kmer = current_kmer;
						ionCorrection2(ii, current_kmer, untrusted_zone_size);
						if(_continue_correction) break;
					}

					
					//end of the read, we should treat this gap here correc snp with trad method
					if(ii == (_readSize-_kmerSize)){
						
						//_bloocoo->__error_detected += 1;
						
						if(!_bloocoo->_ion_mode)
						{
							if(nb_errors_cor == 0){
								nb_checked = _nb_kmers_checked;
								nb_checked =  min(_nb_kmers_checked, untrusted_zone_size);
								nb_checked = max(nb_checked, 0);
								
								nb_errors_cor = aggressiveCorrection(_readSize - untrusted_zone_size - _kmerSize + 1, nb_checked ,RIGHT);
								
								if(nb_errors_cor == 0){
									nb_errors_cor = voteCorrectionInUntrustedZone(_readSize - (untrusted_zone_size-1) - _kmerSize , _readSize-_kmerSize, nb_checked);
								}
								
								if(nb_errors_cor == 0){
									nb_errors_cor = multiMutateVoteCorrectionInUntrustedZone(_readSize - (untrusted_zone_size-1) - _kmerSize, _readSize-_kmerSize, nb_checked,_readSize - untrusted_zone_size ,-1);
								}
							}
						}
					}
					
				}
				
				update_nb_errors_corrected(nb_errors_cor);
				
			} // end of kmers iteration over the read
			
			
		}
		
		if(nb_kmer_trusted == _readSize-_kmerSize+1){
			break;
		}
	}
	
	//if(_bloocoo->_max_trim != 0 && error_exist){
	//	trimSequence();
	//}
	
	writeSequence();
	
	
	if(_bloocoo->_ion_mode   )
	{
//            if(! _newSeq)
//            {
//                _newSeq =  &_readseq ; //new Sequence (_readseq);
//            }
		
//            if(_use_newSeq)
//            {
//                _newSeq->setComment(_sequence->getComment());
//                _newSeq->setIndex(_sequence->getIndex());
//            }
	}
   // printf("fin\n%i %.*s   %i  \n",(int)_newSeq->getDataSize(),(int)_newSeq->getDataSize(),_newSeq->getDataBuffer(),_newSeq->getIndex());

	
	if(_bloocoo->_ion_mode)
	{
	 
		Sequence * outseq;
//            if(_use_newSeq)
//                outseq = _newSeq;
//            else
			outseq =  _sequence ;
		
#ifdef SERIAL
		_outbank->insert (*outseq); //output corrected sequence
#else
		_bankwriter->insert (*outseq);//
		_temp_nb_seq_done ++;
		
		if(_temp_nb_seq_done == 10000) // careful with this val, should be a divisor of buffer size in ordered bank
		{
			_bankwriter->incDone(_temp_nb_seq_done);
			_temp_nb_seq_done=0;
		}
#endif
	   
	}
	else
	{
#ifdef SERIAL
	_outbank->insert (*_sequence); //output corrected sequence
#else
	_bankwriter->insert (*_sequence);
	_temp_nb_seq_done ++;
	
	if(_temp_nb_seq_done == 10000) // careful with this val, should be a divisor of buffer size in ordered bank
	{
		_bankwriter->incDone(_temp_nb_seq_done);
		_temp_nb_seq_done=0;
	}
#endif
	}
	
	
	//    printf("%s\n",s.getDataBuffer());
	//    printf("%s\n",_readseq);
	
	// _bloocoo->_seq_num++; //counter 0 based, not thread-safe  //todo : include the sequence number in the sequence type
	
}
*/



void CorrectReads::execute2(){
	
	int nb_tries =0;
	
	int nb_checked;
	ModelDirect model (_kmerSize);
	ModelDirect::Iterator itKmer (model);
	bool error_exist;
	
	for(int j=0; j<_nb_less_restrictive_correction; j++){
		_nb_kmer_offset += 2;
		
		_continue_correction = true;
		
		while(_continue_correction && nb_tries < 5){
			
			nb_tries += 1; //a remttre à zero avant ce while ? plus lent mais sinon peut shunter les less_restrictive_correction
			
			itKmer.setData (_sequence->getData());
			//itKmer.setData(_sequence->getData());
			//int i = 0;
			
			_kmers.clear();
			for (itKmer.first(); !itKmer.isDone(); itKmer.next()){//, i++){
				_kmers.push_back(itKmer->value());
			}
			
			#ifdef PRINT_DEBUG
			    ModelDirect::Iterator itKmer2 (model);
				kmer_type kmer, kmer_min;
				
				itKmer2.setData (_sequence->getData());
				int i3=0;

				cout << endl << "Read " << _sequence->getIndex() << endl;
				cout << "\t" << _sequence->getDataBuffer() << endl << "\t";
				
				for (itKmer2.first(); !itKmer2.isDone(); itKmer2.next(),i3++){
					kmer = *itKmer2;
					kmer_min = std::min(  revcomp (kmer, _kmerSize),kmer );
					if (_bloom->contains(kmer_min))
						cout << "1";
					else
						cout << "0";
				}
				cout << endl;
			#endif
	
			_continue_correction = searchError(&error_exist);
			//if(correction_finished) break;
		}
		
		if(!error_exist) break;
	}
	
	if(_bloocoo->_max_trim > 0 && error_exist){
		trimSequence();
	}
	

}

void CorrectReads::trimSequence(){
	ModelDirect model (_kmerSize);
	ModelDirect::Iterator itKmer (model);
	kmer_type current_kmer;
	kmer_type current_kmer_min;

	#ifdef PRINT_DEBUG
		ModelDirect::Iterator itKmer2 (model);
		kmer_type kmer, kmer_min;
		
		itKmer2.setData (_sequence->getData());
		int i3=0;

		cout << "\tTrimming Read " << _sequence->getIndex() << endl;
		cout << "\t";
		for (itKmer2.first(); !itKmer2.isDone(); itKmer2.next(),i3++){
			kmer = *itKmer2;
			kmer_min = std::min(  revcomp (kmer, _kmerSize),kmer );
			if (_bloom->contains(kmer_min))
				cout << "1";
			else
				cout << "0";
		}
		cout << endl;
	#endif
	
	int bestNoErrorSize = 0;
	int bestStartPos = 0;
	int noErrorSize = -1;
	int startPos = 0;
	
	itKmer.setData (_sequence->getData());
	
	int i = 0;
	for(itKmer.first(); !itKmer.isDone(); itKmer.next(), i++){
		
		current_kmer = itKmer->value();
		current_kmer_min = min(current_kmer, revcomp(current_kmer, _kmerSize));
		
		if(_bloom->contains(current_kmer_min)){
			if(noErrorSize == -1){
				noErrorSize = 0;
				startPos = i;
			}
			noErrorSize += 1;
		}
		else{
			if(noErrorSize > bestNoErrorSize){
				bestNoErrorSize = noErrorSize;
				bestStartPos = startPos;
			}
			noErrorSize = -1;
		}
	}
	
	if(noErrorSize > bestNoErrorSize){
		bestNoErrorSize = noErrorSize;
		bestStartPos = startPos;
	}
			
	int endPos = bestStartPos + bestNoErrorSize + _kmerSize - 1;
	int trimReadSize = endPos-startPos;
	
	/*
	cout << "best no error size: " << bestNoErrorSize << endl;
	cout << "no error start: " << bestStartPos << endl;
	cout << "no error end: " << endPos << endl;
	cout << "trimmed read size: " << trimReadSize << endl;
	cout << "Can trim: " << (_readSize-trimReadSize <= _bloocoo->_max_trim) << endl;
	*/
	
	memmove(_readseq, _readseq+bestStartPos, trimReadSize);
	_sequence->getData().setSize(trimReadSize);
	//_readseq = trimmed_seq;
}

void CorrectReads::writeSequence(){
	
	/*
	if(_bloocoo->_ion_mode){
	 
		Sequence * outseq;
		//if(_use_newSeq)
		//	outseq = _newSeq;
		//else
		outseq =  _sequence ;
		

	}
	else{
		#ifdef SERIAL
			_outbank->insert (*_sequence); //output corrected sequence
		#else
			_bankwriter->insert (*_sequence);
			_temp_nb_seq_done ++;
			
			if(_temp_nb_seq_done == 10000) // careful with this val, should be a divisor of buffer size in ordered bank
			{
				_bankwriter->incDone(_temp_nb_seq_done);
				_temp_nb_seq_done=0;
			}
		#endif
	}
	*/
	
	#ifdef SERIAL
		_outbank->insert (*_sequence); //output corrected sequence
	#else
		_bankwriter->insert(*_sequence);//
		_temp_nb_seq_done ++;
		
		if(_temp_nb_seq_done == 10000) // careful with this val, should be a divisor of buffer size in ordered bank //was 10000
		{
			_bankwriter->incDone(_temp_nb_seq_done);
			_temp_nb_seq_done=0;
		}
	#endif
}



bool CorrectReads::searchError(bool* error_exist)
{
	*error_exist = false;
	//int _min_successive_trusted_kmers = 2;
	
	//int min_successive_trusted_kmers = _min_successive_trusted_kmers;
	kmer_type current_kmer;
	kmer_type current_kmer_min;
	int nb_errors_cor = 0;
	//int successive_trusted_kmer_count = 0;
				
	#ifdef PRINT_DEBUG
		string debug_str = "\t";
		int debug_str_index = 0;
	#endif
	
	for(int i=0; i<_kmerCount; i++){
	//int i = 0;
	//while(i < _kmerCount){
		
		#ifdef PRINT_DEBUG
			while(debug_str_index < i){
				debug_str += " ";
				debug_str_index += 1;
			}
		#endif
		
		current_kmer = _kmers[i];
		current_kmer_min = min(revcomp(current_kmer, _kmerSize), current_kmer);

		//kmer is solid
		if (_bloom->contains(current_kmer_min)){
			
			#ifdef PRINT_DEBUG
				debug_str += "1";
				debug_str_index += 1;
			#endif
			
			//successive_trusted_kmer_count += 1;
			/*
			if(successive_trusted_kmer_count >= min_successive_trusted_kmers){
				successive_trusted_kmer_count = 0;
				if(i==_kmerCount-1){
					break;
				}
				else if(i==_kmerCount-min_successive_trusted_kmers){
					successive_trusted_kmer_count = 1;
					i += 1;
					continue;
				}
				else{
					i += 0; //_kmerSize/2;//_kmerSize-1;
				}
			}
			
			//Check the last kmers of the read
			if(i>=_kmerCount-1 && successive_trusted_kmer_count < min_successive_trusted_kmers){
				i = _kmerCount-min_successive_trusted_kmers;
				continue;
			}
			*/
		}
		//non solid kmer detected
		else{
			//_bloocoo->__error_detected += 1;
			*error_exist = true;
			//successive_trusted_kmer_count = 0;

			int startPos = searchErrorZoneRec(i, false, 0);
			int endPos = searchErrorZoneRec(i, true, 0);
			
			//int min_successive_possible = _kmerCount - endPos - 1;
			//min_successive_trusted_kmers = min(min_successive_possible, _min_successive_trusted_kmers);
			
			#ifdef PRINT_DEBUG
				debug_str += "0";
				debug_str_index += 1;
				std::string holeDebug = "";
				int holeIndexDebug = 0;
				while(holeIndexDebug<startPos){
					holeDebug += " ";
					holeIndexDebug += 1;
				}
				while(holeIndexDebug<=endPos){
					holeDebug += "0";
					holeIndexDebug += 1;
				}
				printf("\t%s\n", holeDebug.c_str()); 
			#endif
			
			nb_errors_cor += startCorrectionInZone(startPos, endPos);
			i = endPos;
		}
		
		//i += 1;
	}

	#ifdef PRINT_DEBUG
		cout << debug_str << endl;
	#endif
	
	if(nb_errors_cor > 0){
		_local_nb_errors_corrected += nb_errors_cor;
		return true;
	}
	else{
		return false;
	}

}

int CorrectReads::searchErrorZoneRec(int pos, bool extendRight, int trustedKmerCount)
{
	if(extendRight){
		pos += 1;
	}
	else{
		pos -= 1;
	}
	
	if(pos < 0){
		return 0;
	}
	else if(pos > _readSize-_kmerSize){
		return _readSize-_kmerSize;
	}
	
	kmer_type current_kmer = _kmers[pos];
	kmer_type current_kmer_min = min(revcomp(current_kmer, _kmerSize), current_kmer);
		
	if (_bloom->contains(current_kmer_min)){
		trustedKmerCount += 1;
	}
	else{
		trustedKmerCount = 0;
	}
	
	if(trustedKmerCount >= 2){
		if(extendRight){
			return pos-2;
		}
		else{
			return pos+2;
		}
	}
	else{
		return searchErrorZoneRec(pos, extendRight, trustedKmerCount);
	}
	
}

int CorrectReads::startCorrectionInZone(int startPos, int endPos)
{
	int nb_errors_cor = 0;
	int holeSize = endPos - startPos + 1;
	
	#ifdef PRINT_DEBUG
		printf("\t\t[hole_size: %i]  [start_pos: %i]  [end_pos: %i]\n", holeSize, startPos, endPos);
	#endif

	if(holeSize == _kmerSize){
		nb_errors_cor = twoSidedCorrection(endPos);
	}
	
	if(nb_errors_cor > 0) return nb_errors_cor;
	
	int nb_checked = _nb_kmers_checked;
	nb_checked = min(_nb_kmers_checked, holeSize);
	nb_checked = max(nb_checked, 0);	
	
	//Error at begin of read
	if(startPos == 0){
		nb_errors_cor = aggressiveCorrection(endPos, nb_checked, LEFT);
		if(nb_errors_cor > 0) return nb_errors_cor;
		
		nb_errors_cor = voteCorrectionInUntrustedZone(startPos, endPos, nb_checked);
		if(nb_errors_cor > 0) return nb_errors_cor;
		
		nb_errors_cor = multiMutateVoteCorrectionInUntrustedZone(startPos, endPos, nb_checked, -1, endPos);
	}
	//Error at end of read
	else if(endPos == _readSize-_kmerSize){
		nb_errors_cor = aggressiveCorrection(_readSize - holeSize - _kmerSize + 1, nb_checked, RIGHT);
		if(nb_errors_cor > 0) return nb_errors_cor;
		
		nb_errors_cor = voteCorrectionInUntrustedZone(startPos, endPos, nb_checked);
		if(nb_errors_cor > 0) return nb_errors_cor;
		
		//nb_errors_cor = multiMutateVoteCorrectionInUntrustedZone(startPos, endPos, nb_checked, _readSize-holeSize, -1);
		nb_errors_cor = multiMutateVoteCorrectionInUntrustedZone(startPos, endPos, nb_checked, startPos+_kmerSize-1, -1);
								
	}
	//Error in middle of the read
	else{
		nb_errors_cor = aggressiveCorrection(startPos  , nb_checked, RIGHT); // was endPos - _kmerSize + 1
		if(nb_errors_cor > 0) return nb_errors_cor;
		nb_errors_cor = aggressiveCorrection(endPos, nb_checked, LEFT);
		if(nb_errors_cor > 0) return nb_errors_cor;
		
		nb_errors_cor = voteCorrectionInUntrustedZone(startPos, endPos, nb_checked);
		if(nb_errors_cor > 0) return nb_errors_cor;
		
		nb_errors_cor = multiMutateVoteCorrectionInUntrustedZone(startPos, endPos, nb_checked, startPos+_kmerSize-1, endPos);
	}
	

		
	return nb_errors_cor;

}

//not needed if execute2 is used
void CorrectReads::update_nb_errors_corrected(int nb_errors_corrected){
    _local_nb_errors_corrected += nb_errors_corrected;
    if(nb_errors_corrected > 0){
    	_continue_correction = true;
    }
}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
//// twoSidedCorrection
// Method 1 for read correction
// Correct an error in an untrusted zone of size=kmerSize
//pos is the position of the error
int CorrectReads::twoSidedCorrection(int pos)
{
	
	if(!is_pos_correctable(pos)){
		return 0;
	}

	#ifdef PRINT_STATS
		_bloocoo->__correction_methods_calls[Bloocoo::TWO_SIDED] ++;
	#endif
		
    kmer_type left_most_kmer = _kmers[pos-_kmerSize+1];
    kmer_type right_most_kmer = _kmers[pos];
    
    int original_nt = (_readseq[pos]>>1)&3;
    int good_nt;
    kmer_type mutated_kmer, mutated_kmer_min;
    int nb_alternative = 0;
    
	/*
     if(pos > 65 || pos < 33){
     printf("%i\n", pos);
     }*/
	
    for(int nt=0; nt<4; nt++){
		if(nt == original_nt){
			continue;
		}
		
		//Leftmost kmer, we mutate the last NT of this kmer
		mutated_kmer = left_most_kmer;
		mutate_kmer(&mutated_kmer, 0, nt);
		mutated_kmer_min = min(mutated_kmer, revcomp(mutated_kmer, _kmerSize));
		
        
        
		//Check if the mutated left most kmer is trusted or not
        if(_bloom->contains(mutated_kmer_min)){
            
            //The mutated left most kmer is indexed, we can now check if the mutated right most kmer is trusted
            //Rightmost kmer, we mutate the first NT of this kmer
            mutated_kmer = right_most_kmer;
            mutate_kmer(&mutated_kmer, _kmerSize-1, nt);
            mutated_kmer_min = min(mutated_kmer, revcomp(mutated_kmer, _kmerSize));
            
            if(_bloom->contains(mutated_kmer_min)){
            	//Both the left most and right most mutated kmers are trusted so the current NT is a correction alternative
                nb_alternative += 1;
                good_nt = nt;
            }
        }
    }
    
    if(nb_alternative == 1){
    	int nb_correction = apply_correction(pos, good_nt,Bloocoo::TWO_SIDED);
    	#ifdef PRINT_STATS
			_bloocoo->__correction_methods_successes[Bloocoo::TWO_SIDED] += nb_correction;
		#endif
        return nb_correction;
    }
    
    return 0;
    
}



/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
//// aggressiveCorrection
// Method 2 for read correction
// This method can correct an error in an untrusted zone of size > kmerSize (close errors)
//
// nb_kmer_check the number of additional kmers checked
// pos = position of the solid kmer next to the untrusted zone
// kmer_begin  = solid kmer at pos
// Direction:
//     RIGHT: We start from the last trusted kmer BEFORE the untrusted zone and try to correct the first error in this zone.
//     LEFT:  We start from the first trusted kmer AFTER the untrusted zone and try to correct the last error in this zone.
int CorrectReads::aggressiveCorrection(int pos, int nb_kmer_check, Direction direction)
{
    
	//Determine corrected position depending of a left or right agressive correction
    int corrected_pos;
    if(direction ==RIGHT){
        corrected_pos = pos+_kmerSize-1;
    }
    else{
        corrected_pos = pos;
    }
	
	//cout << "\t\t" << _sequence->getIndex() << "    Agressive: [pos: " << corrected_pos << "]" << endl;
    //Cancelled the vote if the corrected position is not correctable
	if(!is_pos_correctable(corrected_pos)){
		return 0;
	}
	
    #ifdef PRINT_STATS
		_bloocoo->__correction_methods_calls[Bloocoo::AGRESSIVE] ++;
	#endif
	
    //Determine the minimum vote threshold
    //It's always the param _nb_min_valid except for the read sides (Sides have limited number of kmers to check)
    int current_max_nb_checkable;
    if(direction ==RIGHT){
    	current_max_nb_checkable = _readSize - corrected_pos;
    }
    else{
    	current_max_nb_checkable = corrected_pos + 1;
    }
    
    int vote_threshold = _nb_min_valid-_nb_kmer_offset;
    if(!_only_decrease_nb_min_valid) nb_kmer_check -= _nb_kmer_offset;
    nb_kmer_check = max(1, nb_kmer_check);
    vote_threshold = max(1, vote_threshold);
    
    //int nb_extra_kmer = max(0, nb_kmer_check-current_max_nb_checkable);
    //nb_extra_kmer = min(1, nb_extra_kmer);
    int nb_extra_kmer = max(0, vote_threshold-current_max_nb_checkable);
    //nb_extra_kmer = min(2, nb_extra_kmer);
    
    //current_max_nb_checkable += nb_extra_kmer;
    //current_max_nb_checkable += nb_extra_kmer;
    //current_max_nb_checkable = nb_extra_kmer + nb_kmer_check;
    
    
    nb_kmer_check = min(current_max_nb_checkable, nb_kmer_check);
    
    
    //If the number of checkable kmers is inferior to the vote threshold then the vote is cancelled
	if(nb_kmer_check+nb_extra_kmer < vote_threshold){
		return 0;
	}
	
	kmer_type kmer_begin = _kmers[pos];
	int votes [4] = {0, 0, 0, 0};
    int good_nt;
    int nb_possibles_firstk=0;
	char nt_temp;
    
    ModelCanonical model (_kmerSize);
    kmer_type current_kmer, current_kmer_min;
    
    int original_nt;
    if(direction ==RIGHT){
        original_nt = (_readseq[corrected_pos]>>1)&3;
    }
    else{
    	original_nt = (_readseq[corrected_pos]>>1)&3;
    }
    
    
	
	
    for(int nt=0; nt<4; nt++)
    {
		if(nt == original_nt){
			continue;
		}
        
		current_kmer = kmer_begin;
		
		if(direction ==RIGHT){
			mutate_kmer(&current_kmer, 0, nt);
		}
		else{
		    mutate_kmer(&current_kmer, _kmerSize-1, nt);
		}
		
		//current_kmer = codeSeedBin(&model, &kmer_begin, nt, direction);
        current_kmer_min = min(current_kmer, revcomp(current_kmer, _kmerSize));
        
        if(_bloom->contains(current_kmer_min)) //kmer is indexed
        {
			votes[nt] += 1;
            nb_possibles_firstk++;
        }
        else{
        	//continue; // if first kmer is not valid, this possible nt is eliminated
        }
        
		//nb_kmer_check-1: -1 because the first mutation above is counted as a kmer check
		for (int ii=0; ii<nb_kmer_check-1; ii++)
		{
		    
		    if(direction == RIGHT){
		        nt_temp = _readseq[corrected_pos+1+ii];
		    }
		    else{
		        nt_temp = _readseq[corrected_pos-1-ii];
		    }
		    
		    codeSeedNT(&model, &current_kmer, nt_temp, direction);
		    
		    current_kmer_min = min(current_kmer, revcomp(current_kmer, _kmerSize));
		    if(_bloom->contains(current_kmer_min)) //kmer is indexed
		    {
		        votes[nt] += 1;
		    }
		    
		}
		
		if(nb_extra_kmer > 0 /*&& votes[nt] > 0*/){
			extendedAgressiveCorrection(votes, &model, &current_kmer, nt, nb_extra_kmer, direction);
		}
    }
    
    
    //print_agressive_votes(votes);
    
	//Searching max score in votes
	int max_score = votes[0];
	for(int i=1; i<4; i++){
 		if(votes[i] > max_score){
			max_score = votes[i];
		}
	}
    
	//int t = (nb_kmer_check*25)/100;
	//int t = max(1, nb_kmer_check);
	//int t = 1; //max(1, nb_kmer_check);
	//if(nb_possibles_firstk==1) vote_threshold=min(current_max_nb_checkable,3);
    
	if(max_score < vote_threshold){
		return 0;
	}
	
	//printf("max: %i\n", max_score);
    
	//We have the max score in votes, now we have to check if this score
	//is uniq or if more than one nt has the max score
	int nb_alternative = 0;
    
	for(int i=0; i<4; i++){
 		if(votes[i] == max_score){
			nb_alternative += 1;
			good_nt = i;
		}
	}
    
	//printf("good_nt: %i, nb_alternative: %i\n", good_nt, nb_alternative);
    
	//if max score is uniq we can correct the sequence
	//else if nb_nt_max > 1 then there are more than one candidate for the correction
	//so we cannot determine a good correction, the correction is cancelled.
	if(nb_alternative == 1){
		int nb_correction = apply_correction(corrected_pos, good_nt, Bloocoo::AGRESSIVE);
    	#ifdef PRINT_STATS
			_bloocoo->__correction_methods_successes[Bloocoo::AGRESSIVE] += nb_correction;
		#endif
        return nb_correction;
	}
    
    return 0;
}



/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
int CorrectReads::voteCorrectionInUntrustedZone(int start_pos, int end_pos, int nb_kmer_checked){
	
    #ifdef PRINT_STATS
		_bloocoo->__correction_methods_calls[Bloocoo::VOTE] ++;
	#endif
	
	//printf("%i:    %i %i\n", _seq_num, start_pos, end_pos);
	//printf("\n\t");
	//(*kmers[start_pos]).printASCII(_kmerSize);
	//(*kmers[end_pos]).printASCII(_kmerSize);
	
	//start_pos = max(0, start_pos);
	
	int nb_errors_cor = 0;
    
	while(nb_errors_cor == 0 && start_pos < end_pos){
		int untrusted_zone_size = end_pos - start_pos;
		int new_nb_checked = min(nb_kmer_checked, untrusted_zone_size);
		if(new_nb_checked <= 0){ break; }
		
		nb_errors_cor = voteCorrection(start_pos, new_nb_checked);
		start_pos += _kmerSize/2;//new_nb_checked;
	}
	
	return nb_errors_cor;
}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
//// voteCorrection
// Method 3 for read correction
// This method is used if the method 2 fails to correct an error
//
// nb_kmer_check the number of additional kmers checked
// pos: the last trusted position in readseq before the untrusted zone
// kmer_begin: solid kmer at pos, it's the last trusted kmer before the untrusted zone
int CorrectReads::voteCorrection(int start_pos, int nb_kmer_check)
{

	
	int vote_threshold = _nb_min_valid - _nb_kmer_offset;
	if(!_only_decrease_nb_min_valid) nb_kmer_check -= _nb_kmer_offset;
	nb_kmer_check = max(1, nb_kmer_check);
	vote_threshold = max(1, vote_threshold);
	
	if(nb_kmer_check < vote_threshold){
		return 0;
	}
	
    int good_nt;
    //KmerModel model (_kmerSize);
    //kmer_type kmer_begin = *kmers[start_pos];
    kmer_type current_kmer, current_kmer_min, mutated_kmer;
    
	int nb_column = nb_kmer_check + _kmerSize - 1;
	int votes[nb_column][4];
	for (int i = 0; i < nb_column; i++) {
		for (int nt=0; nt < 4; nt++) {
			votes[i][nt] = 0;
		}
	}
	
	//kmer_begin.printASCII(_kmerSize);
	//current_kmer = kmer_begin;
	
	for(int i=0; i < nb_kmer_check; i++){
        
		int read_pos = start_pos + i;
	
		current_kmer = _kmers[read_pos];
		//current_kmer = model.codeSeedRight(current_kmer, _readseq[read_pos+(_kmerSize-1)], Data::ASCII, KMER_DIRECT);
		
		//printf("\tcurrent kmer:    "); current_kmer.printASCII(_kmerSize);
		
		current_kmer_min = min(current_kmer, revcomp(current_kmer, _kmerSize));
		//if(_bloom->contains(current_kmer_min)){
		//    continue;
		//}
		//if(!is_pos_correctable(read_pos, _readseq)){
		//	continue;
		//}
        
		for(int kpos=0; kpos < _kmerSize; kpos++){
			if(!is_pos_correctable(read_pos+kpos)){
				continue;
			}
            
			int original_nt = (_readseq[read_pos+kpos]>>1)&3;
			
    		for(int nt=0; nt<4; nt++){
				if(nt == original_nt){
					continue;
				}
				
				mutated_kmer = current_kmer;
				mutate_kmer(&mutated_kmer, _kmerSize-1-kpos, nt);
				//mutated_kmer.printASCII(_kmerSize);
				
				/*
                 printf("mutate at %i with %c\n", _kmerSize-1-kpos, bin2NT[nt]);
                 printf("curre_kmer: ");
                 current_kmer.printASCII(_kmerSize);
                 printf("mutated_kmer: ");
                 mutated_kmer.printASCII(_kmerSize);
                 printf("\n");
                 */
				
				current_kmer_min = min(mutated_kmer, revcomp(mutated_kmer, _kmerSize));
				if(_bloom->contains(current_kmer_min))
				{
				    votes[i+kpos][nt] += 1;
					//print_votes(votes);
				}
    		}
		}
	}
    
	//print_votes(votes, nb_column);
    
	//Search max vote in the matrix
	int vote, maxVote = 0 , nb_max_vote= 0;
	for (int i = 0; i < nb_column; i++) {
		for (int nt=0; nt < 4; nt++) {
			vote = votes[i][nt];
            if (vote == maxVote) {
                nb_max_vote ++ ;
			}
			if (vote > maxVote) {
				maxVote = vote;
                nb_max_vote =1;
			}
		}
	}
    
    
    if(nb_max_vote !=1 ){
		return 0;
	}
    
	if(maxVote < vote_threshold){
		return 0;
	}
	
	int nb_alternative, corrected_pos;
	int nb_correction = 0;
    
	
    //si plusieurs max vote, on les corrige tous : ptet cest pas bien : ils sont peut etre plus proches qun kmer, normalement impossible avec  1 muta
	//For each column in the matrix, if the max vote is find and is uniq we can apply the correction
	for (int i = 0; i < nb_column; i++) {
		nb_alternative = 0;
		for (int nt=0; nt < 4; nt++) {
			vote = votes[i][nt];
			if (vote == maxVote) {
				nb_alternative += 1;
				good_nt = nt;
			}
		}
		if(nb_alternative == 1){
			corrected_pos = start_pos + i;
			int nb_cor = apply_correction(corrected_pos, good_nt,Bloocoo::VOTE);
			nb_correction += nb_cor;
	    	#ifdef PRINT_STATS
				_bloocoo->__correction_methods_successes[Bloocoo::VOTE] += nb_cor;
			#endif
		}
	}
	
    return nb_correction;
}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
int CorrectReads::apply_correction(int pos, int good_nt,int algo){
	if(!is_pos_correctable(pos)){
		return 0;
	}
	
    static const char* methodname[] = { "twosided", "aggressive", "vote", "multi" , "side"  };
    

	if(_wantErrTabFile){
	    LocalSynchronizer synchro (_bloocoo->getErrtabMutex());
	    _bloocoo->_errfile->print("%i\t%i\t%c\t%c\n",_sequence->getIndex(),pos,bin2NT[good_nt],_readseq[pos]);
	    _bloocoo->_errfile_full->print("%i\t%i\t%c\t%c;%s\n",_sequence->getIndex(),pos,bin2NT[good_nt],_readseq[pos],methodname[algo]);
    }
    
    _readseq[pos] = bin2NT[good_nt]; //correc sequence in ram
    _corrected_pos.push_back(pos);
    return 1;
}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
void CorrectReads::codeSeedBin(ModelCanonical* model, kmer_type* kmer, int nt, Direction direction){
	
	if(direction == RIGHT)
	{
		/** We initialize the canonical kmer. */
		ModelCanonical::Kmer tmp;  tmp.set (*kmer, revcomp(*kmer, _kmerSize));
		
		/** We get the kmer successor (in particular its value). */
		tmp = model->codeSeedRight (tmp, nt, Data::INTEGER); //.value()
		
		*kmer = tmp.forward() ;  // GR I think we need tmp.forward(), was .value before
		
	}
	else
	{
		/** We initialize the canonical kmer. */
		ModelCanonical::Kmer tmp;  tmp.set (revcomp(*kmer, _kmerSize), *kmer);
		
		/** We get the kmer successor. */
		tmp = model->codeSeedRight (tmp, binrev[nt], Data::INTEGER);
		
		*kmer = tmp.forward() ;
		/** We get the wanted value. */
		//  *kmer = tmp.value()==tmp.forward() ? tmp.revcomp() : tmp.forward(); // GR hmm no I think we need tmp.forward()
	}
}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
void CorrectReads::codeSeedNT(ModelCanonical* model, kmer_type* kmer, char nt, Direction direction){
	//Attention ici, penser à utiliser une fonction binToNT statique.
	return codeSeedBin(model, kmer, (nt>>1)&3, direction);
}


/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
//// mutate_kmer
// fonction to mutate kmer : takes kmer, pos (de la fin ,0-based), and nt (nt = 0,1,2ou 3)
// par exemple  kmer ,2 , C    avec kmer =  AAAAAAAAAA
//
// return :
// AAAAAAACAA
void CorrectReads::mutate_kmer(kmer_type * kmer, int pos, char nt)
{
    kmer_type trois; trois.setVal(3);
    kmer_type reset_mask =  ~(trois << (pos*2));
    kmer_type kmer_nt; kmer_nt.setVal(nt);
    kmer_type set_mask =  kmer_nt << (pos*2);
    *kmer = (*kmer & reset_mask )  | set_mask;
}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
bool CorrectReads::is_pos_correctable(int pos){
	bool N_at_pos = (_readseq[pos] == 'N');
	bool is_pos_corrected = (std::find(_corrected_pos.begin(), _corrected_pos.end(), pos) != _corrected_pos.end());
	return !N_at_pos && !is_pos_corrected;
}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
int CorrectReads::multiMutateVoteCorrectionInUntrustedZone(int start_pos, int end_pos, int nb_kmer_checked, int expected_first_pos, int expexted_second_pos){
	
    #ifdef PRINT_STATS
		_bloocoo->__correction_methods_calls[Bloocoo::MULTI_MUTATE_VOTE] ++;
	#endif
	
	//printf("%i:    %i %i\n", _seq_num, start_pos, end_pos);
	//printf("\n\t");
	//(*kmers[start_pos]).printASCII(_kmerSize);
	//(*kmers[end_pos]).printASCII(_kmerSize);
	
	//start_pos = max(0, start_pos);
	//printf("\n\tMulti mutate start !!!!!!!!!!\n");
	
	//_tab_multivote = (unsigned char *) malloc(TAB_MULTIVOTE_SIZE*sizeof(unsigned char)); // pair of muta  = 16 nt *128 pos * 16 (max dist)

    
	#if PRINT_LOG_MULTI
		_debug->print("%i\t",end_pos-start_pos +1);
	#endif
	int nb_errors_cor = 0;
	int nz=0;
	while(nb_errors_cor == 0 && start_pos < end_pos){
		int untrusted_zone_size = end_pos - start_pos;
		int new_nb_checked = min(nb_kmer_checked, untrusted_zone_size);
		if(new_nb_checked <= 0){ break; }
		
        nz++;
		nb_errors_cor = multiMutateVoteCorrection(start_pos, new_nb_checked,expected_first_pos,expexted_second_pos);
		start_pos += new_nb_checked;
		// start_pos += _kmerSize/2;
	}
	#if PRINT_LOG_MULTI
		_debug->print("%i\t%i\n",nz,nb_errors_cor);
	#endif
	return nb_errors_cor;
}


/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
int CorrectReads::multiMutateVoteCorrection(int start_pos, int nb_kmer_check,int expected_first_pos, int expected_second_pos)
{
	
    int current_max_nb_checkable = 999999;
    if(expected_second_pos > 0 &&  expected_second_pos < 10){
        current_max_nb_checkable = 3;
    }
    if(expected_first_pos > 0   &&  (_readSize -expected_first_pos)<10 ){
    	current_max_nb_checkable = 3;
    }
    
    int vote_threshold = min(current_max_nb_checkable, _nb_min_valid - _nb_kmer_offset);
    
	//int vote_threshold = _nb_min_valid-kmer_offset;
	if(!_only_decrease_nb_min_valid) nb_kmer_check -= _nb_kmer_offset;
	nb_kmer_check = max(1, nb_kmer_check);
	vote_threshold = max(1, vote_threshold);
	
    
	if(nb_kmer_check < vote_threshold){
		return 0;
	}
	
	memset(_tab_multivote, 0, TAB_MULTIVOTE_SIZE*sizeof(unsigned char));
    
	int nb_max_vote = 0;
	int max_vote = 0;
	int vote;
	int nb_correction = 0;
	int good_index;
	
	//int nb_column = nb_kmer_check + _kmerSize - 1;
	//std::map<std::string, int> votes;
	
	for(int i=0; i < nb_kmer_check; i++){
		kmer_type current_kmer = _kmers[start_pos+i];
		multiMutateVoteCorrectionRec(start_pos, i, current_kmer, nb_kmer_check, 0, 0, &max_vote, &nb_max_vote, &good_index, 0, 0,expected_first_pos,expected_second_pos);
		//multiMutateVoteCorrectionRec(start_pos, i, current_kmer, kmers, nb_kmer_check, 0, 0, _tab_multivote, 0, 0);
	}
    
	
	if(max_vote < vote_threshold){
		return 0;
	}
	
	
	if(nb_max_vote != 1){
		return 0;
	}
	
    
	int pos1, dist, nt1, nt2;
	decode_index(good_index, &pos1, &dist, &nt1, &nt2);
	
	int nb_cor = apply_correction(start_pos+pos1, nt1,Bloocoo::MULTI_MUTATE_VOTE);
	#ifdef PRINT_STATS
		_bloocoo->__correction_methods_successes[Bloocoo::MULTI_MUTATE_VOTE] += nb_cor;
	#endif
	nb_correction += nb_cor;
    
	nb_cor = apply_correction(start_pos+pos1+dist, nt2,Bloocoo::MULTI_MUTATE_VOTE);
	#ifdef PRINT_STATS
		_bloocoo->__correction_methods_successes[Bloocoo::MULTI_MUTATE_VOTE] += nb_cor;
	#endif
	nb_correction += nb_cor;
    
	return nb_correction;
    
}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
void CorrectReads::multiMutateVoteCorrectionRec(int start_pos, int kmer_offset, kmer_type current_kmer, int nb_kmer_check, int kmer_index, int current_nb_mutation, int* max_vote, int* nb_max_vote, int *good_index, int pos1, int nt1, int expected_first_pos, int expected_second_pos){
	int read_pos = start_pos + kmer_offset;
	
	
	//Determine the distance between 2 mutations
	int end_kpos = _kmerSize;
	if(current_nb_mutation == 1){
		end_kpos = std::min(kmer_index+_max_multimutation_distance, (int)_kmerSize);
		//printf("%i %i\n", kmer_index, end_kpos);
	}
	
	//Iterating on kmer positions from kmer_index to the max mutate distance
	//If the max number of mutation is not reached, we iterate over all the kmer positions
	for(int kpos=kmer_index; kpos < end_kpos; kpos++){
		
		//Don't vote if the current position is not correctable
		if(!is_pos_correctable(read_pos+kpos)){
			continue;
		}
        
        //si celui ci a abs 1
        //expected first pos is known, restrict search to +-1
		if((expected_first_pos >=0) &&
           (current_nb_mutation == 0) &&
           abs(read_pos+kpos - expected_first_pos) > 0 ){ //only test pos +-1 of expected error pos
			continue;
		}
        
        //ceux là (extension vers left) ne semblent pas impacter recall, et gain perf correct
        //expected second pos is known  restrict first to max_dist away
        if((expected_second_pos >=0) &&
           (current_nb_mutation == 0) &&
           abs(read_pos+kpos - expected_second_pos) >  (1+_max_multimutation_distance) ){
			continue;
		}
        
        if((expected_second_pos >=0) &&
           (current_nb_mutation == 1) &&
           abs(read_pos+kpos - expected_second_pos) > 0 ){
			continue;
		}
        
        
		int original_nt = (_readseq[read_pos+kpos]>>1)&3;
		
		//Iterating A, C, G, T except the original_nt
		for(int nt=0; nt<4; nt++){
			if(nt == original_nt){
				continue;
			}
			
			//Mutate the current kmer
			kmer_type mutated_kmer = current_kmer;
			mutate_kmer(&mutated_kmer, _kmerSize-1-kpos, nt);
			
			/*
             //Add this new mutation ("pos nt") to the string that stores the mutations chain. (Ex: "60 A 12 B...")
             //_oss.clear();
             _oss.str("");
             _oss << kmer_offset+kpos;
             std::string mutations_copy = mutations + " " + _oss.str();
             //_oss.clear();
             _oss.str("");
             _oss << nt;
             mutations_copy += " " + _oss.str();
             */
			
			//If the max number of mutations is reached, we can start voting for the mutated kmers
			if(current_nb_mutation == 1){
				/*
                 printf("kmer index: %i\n", kpos);
                 printf("\tcurent kmer:  "); current_kmer.printASCII(_kmerSize);
                 printf("\tmutated kmer:  "); mutated_kmer.printASCII(_kmerSize);
                 printf("\n");*/
				
				kmer_type mutated_kmer_min = min(mutated_kmer, revcomp(mutated_kmer, _kmerSize));
				if(_bloom->contains(mutated_kmer_min)){
					//int dist = kpos - kmer_index +1;
					//printf("%i %c %i %c\n", pos1+start_pos, bin2NT[nt1], dist, bin2NT[nt]);
					int index = make_index(pos1, kpos-kmer_index+1, nt1, nt);
					_tab_multivote[index] += 1;
					
					if(_tab_multivote[index] == *max_vote){
						*nb_max_vote += 1;
					}
					else if(_tab_multivote[index] > *max_vote){
						*max_vote = _tab_multivote[index];
						*nb_max_vote = 1;
						*good_index = index;
					}
					
                    
				}
			}
			//If the max number of mutations is not reached, we have to recurcivelly continue to applied mutations before voting
			//The kmer_index and the current number of mutation is incremented 
			else{
				multiMutateVoteCorrectionRec(start_pos, kmer_offset, mutated_kmer, nb_kmer_check, kpos+kmer_index+1, current_nb_mutation+1, max_vote, nb_max_vote, good_index, kmer_offset+kpos, nt,expected_first_pos,expected_second_pos);
			}
		}
	}
}




/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
void CorrectReads::extendedAgressiveCorrection(int votes[4], ModelCanonical* model, kmer_type* last_mutated_kmer, int mutated_nt, int nb_kmer_check, Direction direction)
{
    //#ifdef PRINT_STATS
	//	_bloocoo->__correction_methods_calls[Bloocoo::READ_SIDE] ++;
	//#endif
    
	bool posVoted[nb_kmer_check];
	memset(posVoted, 0, nb_kmer_check);
	
	//for(int i=0; i<nb_kmer_check; i++){
	//	posVoted[i] = false;
	//}
	
	//printf("%i\n", posHasVote[0]);
	//for(int i=0; i<nb_kmer_check; i++){
    
	//printf("------------------------------------------------\n");
	//printf("last kmer: "); (*last_mutated_kmer).printASCII(_kmerSize);
	//print_agressive_votes(votes);
	extendedAgressiveCorrectionRec(votes, model, last_mutated_kmer, mutated_nt, nb_kmer_check, direction, posVoted, 0);
	//print_agressive_votes(votes);
	//printf("------------------------------------------------\n");
}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
void CorrectReads::extendedAgressiveCorrectionRec(int votes[4], ModelCanonical* model, kmer_type* mutated_kmer, int mutated_nt, int nb_kmer_check, Direction direction, bool posVoted[], int depth)
{
	if(depth >= nb_kmer_check){
		return;
	}
	
	for(int nt=0; nt<4; nt++){
		
		kmer_type current_kmer = *mutated_kmer;
		codeSeedBin(model, &current_kmer, nt, direction);
		
		
		
		kmer_type current_kmer_min = min(current_kmer, revcomp(current_kmer, _kmerSize));
		if(_bloom->contains(current_kmer_min)) //kmer is indexed
		{
			if(!posVoted[depth]){
				//current_kmer.printASCII(_kmerSize);
				votes[mutated_nt] += 1;
				posVoted[depth] = true;
			}
			
			extendedAgressiveCorrectionRec(votes, model, &current_kmer, mutated_nt, nb_kmer_check, direction, posVoted, depth+1);
			
		}
	}
    
}



/*********************************************************************
 ** METHOD  :make_indice
 ** PURPOSE :encode mutation pair info (pos1, dist_to_pos2, nt1,nt2) into a single int
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
unsigned int CorrectReads::make_index(int pos1, int dist, int nt1, int nt2)
{
    unsigned int idx = pos1;
    
    idx <<= 4;
    idx |=  (dist & 15); //dist is constrained to [0-15]
    idx <<= 2;
    idx |=  nt1; 
    idx <<= 2;
    idx |=  nt2;
    
    return idx;
}


/*********************************************************************
 ** METHOD  :make_indice
 ** PURPOSE :encode mutation pair info (pos1, dist_to_pos2, nt1,nt2) into a single int
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
int CorrectReads::decode_index(unsigned int idx, int * pos1, int * dist, int * nt1, int * nt2)
{
    *nt2 = idx & 3;
    idx >>= 2;
    *nt1 = idx & 3;
    idx >>= 2;
    *dist = idx & 15;
    idx >>= 4;
    *pos1 = idx;
    
    return *pos1;
}




/*
//// print_agressive_votes
//Debug function for agressive correction, print the votes tab
void Bloocoo::print_agressive_votes(int votes[4]){
	printf("#############################################\n");
	for(int i=0; i<4; i++){
		printf("%c: %i\n", bin2NT[i], votes[i]);
	}
	printf("#############################################\n");
}

//print_votes
//Debug function for vote correction, print the votes matrix
void Bloocoo::print_votes(int votes[][4], int nb_column){
	printf("#############################################\n");
	for(int nt=0; nt<4; nt++){
		printf("%c  ", bin2NT[nt]);
		for(int i=0; i<nb_column; i++){
			printf("%i ", votes[i][nt]);
		}
		printf("\n");
	}
	printf("#############################################\n");
}
void Bloocoo::print_read_correction_state(KmerModel* model, Sequence& s){
	KmerModel::Iterator itKmer2 (*model);
	kmer_type kmer, kmer_min;
	
	itKmer2.setData (s.getData());
	int i3=0;
	//printf("\nread correction state (%i):\n", seq_num);
	//printf("\t%s\n", s.getDataBuffer());
	//printf("\t");
	
	std::ostringstream oss;
    oss << s.getIndex();
    
    
	__badReadStack += "Read correction state (" + oss.str() += "):\n";
	__badReadStack += "\t" + std::string(s.getDataBuffer());
	__badReadStack += "\n\t";
	
	for (itKmer2.first(); !itKmer2.isDone(); itKmer2.next(),i3++){
		kmer = *itKmer2;
		kmer_min = std::min(  revcomp (kmer, _kmerSize),kmer );
		if (_bloom->contains(kmer_min)){
			//printf("1");
			__badReadStack += "1";
		}
		else{
			//printf("0");
			__badReadStack += "0";
		}
	}
	//printf("\n");
	__badReadStack += "\n";
}

void Bloocoo::print_read_if_not_fully_corrected(KmerModel* model, Sequence& s){
	std::string read_state = "";
	std::string normal_state = "";
	KmerModel::Iterator itKmer2 (*model);
	kmer_type kmer, kmer_min;
	
	itKmer2.setData (s.getData());
	int i3=0;
	//printf("\nread correction state (%i):\n", _seq_num);
	//printf("\t%s\n", s.getDataBuffer());
	//printf("\t");
	for (itKmer2.first(); !itKmer2.isDone(); itKmer2.next(),i3++){
		normal_state += "1";
		kmer = *itKmer2;
		kmer_min = std::min(  revcomp (kmer, _kmerSize),kmer );
		if (_bloom->contains(kmer_min)){
			read_state += "1";
		}
		else{
			read_state += "0";
		}
	}
	if(normal_state.compare(read_state) != 0){
		//print_read_correction_state(_bloom, model, _kmerSize, s, _seq_num);
		printf("%s", __badReadStack.c_str());
	}
	
}
*/



void CorrectReads::executeIonCorrection(){
	
	int nb_tries =0;
	
	
	
	int nb_checked;
	ModelDirect model (_kmerSize);
	ModelDirect::Iterator itKmer (model);
	kmer_type current_kmer;
	kmer_type current_kmer_min;
	kmer_type last_wrong_kmer; last_wrong_kmer.setVal(0); // else it is not initialized
	int pos_homopo = 0;
	

	
	for(int j=0; j<_nb_less_restrictive_correction; j++){
		_nb_kmer_offset += 2;
		
		
		int nb_kmer_trusted = 0;
		_continue_correction = true;
		bool first_gap = true;
		int ii=0;
		
		
		while(_continue_correction && nb_tries <5)
		{
			nb_tries ++;
			nb_kmer_trusted = 0;
			_continue_correction = false;
			first_gap = true;
			ii = 0;
			
			itKmer.setData (_sequence->getData());
			_readSize = _sequence->getDataSize();
			
#ifdef DEBION
			 printf("%.*s   %i  \n",(int)_sequence->getDataSize(),_sequence->getDataBuffer(),_sequence->getIndex());
#endif
			int untrusted_zone_size = 0;
			int trusted_zone_size = 0;
			int real_untrusted_zone_size = 0;
			
			//Mettre en dehors du while dans une version final (attention dangereux), faire gaffe a ne jamais utiliser un kmer de ce tableau
			//avec un indice > ii
			//int ks = _readSize-_kmerSize+1;
			//if (ks <= 0) ks =10;
			//kmer_type* kmers[ks]; // there was a bug with [_readSize-_kmerSize+1], array should always be init with >0 size
			_kmers.clear();
			
			
			//if(PRINT_DEBUG){ _bloocoo->print_read_correction_state(&model, _sequence);}
			
			// getSynchro()->lock()  ;
			pos_homopo=0;
			
			// We iterate the kmers of this sequence
			for (itKmer.first(); !itKmer.isDone(); itKmer.next(),ii++)
			{
				
				
				
				
				int nb_errors_cor = 0;
				
				current_kmer = itKmer->value();
				//kmers[ii] = &(*itKmer);
				_kmers.push_back(current_kmer);
				
				current_kmer_min = min(revcomp(current_kmer, _kmerSize), current_kmer);
				
				
				
				//Si on veut checker le dernier kmer du read et qu'on est dans une zone untrusted alors on considere
				//automatiquement ce dernier kmer comme untrusted et on laisse les démarches de correction de fin de read s'effectuer.
				bool is_last_kmer_indexed_after_hole = (ii==_readSize-_kmerSize && trusted_zone_size==0);
				
				
				if (!is_last_kmer_indexed_after_hole && _bloom->contains(current_kmer_min)) //kmer is solid
				{
					
					nb_kmer_trusted += 1;
					
					trusted_zone_size += 1;
					
					//beginning of indexed zone
					if(trusted_zone_size == 2)
					{
						//printf("pos 2 pos_homopo %i \n",pos_homopo);
						//todo put this in outside func
						int insert_pos =ii-2;
						int dele_pos= ii-1;

						if(pos_homopo==1 && insert_pos >=5) // todo :le prob des bouts
						{
							
							kmer_type corrected_kmer = last_wrong_kmer;
							//il faut faire confirmation avec kmer solid qd meme ... sinon boucle insert/dele et rapetisse read
							#ifdef DEBION
							printf("insert M %i lenm %i rl %i\n",ii-2,_readSize-(ii-2)-1,_readSize);
#endif
						   // std::cout << " last wrong kmer    " << last_wrong_kmer.toString(_kmerSize) << std::endl;
							//il faut rtempalcer premiere nt par celle davant
							
							int ins_size = 1;

							for(ins_size = 1; ins_size<=3;ins_size ++ )
							{
								corrected_kmer = last_wrong_kmer;
								mutate_kmer(&corrected_kmer, _kmerSize-1, NT2int(_readseq[ii-2-ins_size]));
							 //   if((ii-2-ins_size)>=_sequence->getDataSize()  || (ii-2-ins_size) < 0 )
							  //      printf("PB %i %i %i get %i   kk %i \n",insert_pos,ins_size,_readSize,_sequence->getDataSize(),ii-2-ins_size);
								
							//    std::cout << "try corrected kmer    " << corrected_kmer.toString(_kmerSize)  <<" ins size" << ins_size << std::endl;

								if (_bloom->contains(min(revcomp(corrected_kmer, _kmerSize), corrected_kmer)))
								{

									memmove(_readseq+insert_pos + 1 - ins_size , _readseq + insert_pos + 1, _readSize-(insert_pos)-1);
									_continue_correction = true;
									_sequence->getData().setSize(_readSize-ins_size);
									_local_nb_ins_corrected += ins_size;
#ifdef DEBION
									printf("mm %i %i  %i\n",+insert_pos + 1 - ins_size ,  + insert_pos + 1, _readSize-(insert_pos)-1);

									printf("apply correc insert len %i pos %i \n",ins_size,insert_pos);
#endif

									break;
									
								}
							}


						}
						else if (pos_homopo==2 && dele_pos >=5) // todo :le prob des bouts
						{
							kmer_type corrected_kmer = last_wrong_kmer;

						 //   printf("possible deletion in read  at pos %i  len %i  _readSize %i \n",ii-1,_readSize-(ii-1)-1,_readSize);
#ifdef DEBION

							printf("dele M %i lenm %i rl %i\n",ii-1,_readSize-(ii-1)-1,_readSize);
						  //  printf("%c %c \n",_readseq[dele_pos],_readseq[dele_pos-1]);
						  //  std::cout << " last wrong kmer    " << last_wrong_kmer.toString(_kmerSize) << std::endl;
#endif

							
							int del_size = 1;
							
							for(del_size = 1; del_size<=3;del_size ++ )
							{
								corrected_kmer = last_wrong_kmer;
								corrected_kmer = corrected_kmer >> (2*del_size);
								mutate_kmer(&corrected_kmer, _kmerSize-1, NT2int(_readseq[dele_pos-1]));
								for(int hh=0; hh<del_size; hh++)
								{
									mutate_kmer(&corrected_kmer, _kmerSize-2-hh, NT2int(_readseq[dele_pos]));
								}
							 //   std::cout << "try corrected kmer    " << corrected_kmer.toString(_kmerSize) <<" del size" << del_size<<   std::endl;
								
								
								if (_bloom->contains(min(revcomp(corrected_kmer, _kmerSize), corrected_kmer)))
								{
									memmove(_readseq+dele_pos+del_size,_readseq+dele_pos,_readSize-dele_pos-del_size);
									_continue_correction = true;
									_sequence->getData().setSize(_readSize);
									_local_nb_del_corrected += del_size;
#ifdef DEBION
									printf("mm %i %i  %i\n",dele_pos+del_size ,  + dele_pos, _readSize-dele_pos-del_size);

									printf("apply correc del len %i pos %i \n",del_size,dele_pos);
#endif
									break;
									
								}
								
							}

						}
						if (_continue_correction) break;
						
					}
					
					if(trusted_zone_size==2){
						first_gap = false;
					}
					
					////////////traitement d'un faux positif du bloom isolé
					//   if(trusted_zone_size==1) tai_previous_break = untrusted_zone_size; // begin of indexed zone    ( if untrusted_zone_size ?)
					if (trusted_zone_size > 1) // do not reset  untrusted_zone_size if a single positive (might be a FP)
					{
						untrusted_zone_size = 0; //reset previous hole size
					}
					if (trusted_zone_size==1)// begin of indexed zone
					{
						untrusted_zone_size++; //  treat a single solidkmer  as an erroneous kmer
					}
					//////////
					
					//version sans traitement du FP du bloom
					//if(trusted_zone_size==1){
					//	untrusted_zone_size = 0;
					//}
					
				}
				else //kmer contains an error
				{
					untrusted_zone_size += 1;
					trusted_zone_size = 0;
					
					last_wrong_kmer = current_kmer;
					
					//todo put this in outside func
					pos_homopo=contains_homopolymer(current_kmer,_kmerSize);
					//
//                            if (pos_homopo)
//                            {
//                                std::cout << "    " << current_kmer.toString(_kmerSize)<< "  " << pos_homopo << " " << ii<< std::endl;
//                            }

					
					if(untrusted_zone_size==1)//first wrong kmer
					{
						pos_homopo=contains_homopolymer( revcomp(current_kmer, _kmerSize),_kmerSize);
						
						kmer_type wrong_kmer = revcomp(current_kmer, _kmerSize);

					 //   std::cout << "  wrong kmer    " << wrong_kmer.toString(_kmerSize) << std::endl;

						kmer_type corrected_kmer;
						
						//printf("frist wrong kmer homopo %i pos %i\n",pos_homopo,ii+_kmerSize-1);
						int insert_pos = ii+_kmerSize-1;
						int delete_pos = ii+_kmerSize-2;

						if(pos_homopo==1 &&  insert_pos >  (_readSize-_kmerSize) ) //&& insert_pos < _readSize-5
						{
						   //   printf("insert F  %i lenm %i rl %i\n",insert_pos,_readSize-(insert_pos)-1,_readSize);

							
							int ins_size = 1;
							
							for(ins_size = 1; ins_size<=3;ins_size ++ )
							{
								corrected_kmer = wrong_kmer;
								mutate_kmer(&corrected_kmer, _kmerSize-1, binrev[NT2int(_readseq[insert_pos+ins_size])]);
							//     std::cout << " corrected kmer    " << corrected_kmer.toString(_kmerSize) << std::endl;
								
								if (_bloom->contains(min(revcomp(corrected_kmer, _kmerSize), corrected_kmer)))
								{
								memmove(_readseq+insert_pos+1-ins_size, _readseq+insert_pos+1, _readSize-(insert_pos)-1);
								_continue_correction = true;
								_sequence->getData().setSize(_readSize-ins_size);
								_local_nb_ins_corrected += ins_size;
								//    printf("apply correc ins F  len %i pos %i\n",ins_size,insert_pos);

								break;
								}
								
							}

						}
						
						if(pos_homopo==2 &&  delete_pos >  (_readSize-_kmerSize) && delete_pos < _readSize-5) //ne gere pas la toute fin // todo :le prob des bouts
						{

#ifdef DEBION
							printf("dele  F  %i lenm %i  rl %i\n",delete_pos,_readSize-delete_pos-1,_readSize);
#endif
							
							int del_size = 1;

							for(del_size = 1; del_size<=3;del_size ++ )
							{
								corrected_kmer = wrong_kmer;
								//std::cout << " corrected kmer  oo  " << corrected_kmer.toString(_kmerSize)  <<  std::endl;

								//create the correcetd kmer
								corrected_kmer = corrected_kmer >> (2*del_size);
								mutate_kmer(&corrected_kmer, _kmerSize-1, binrev[NT2int(_readseq[delete_pos+1])]);
								for(int hh=0; hh<del_size; hh++)
								{
									mutate_kmer(&corrected_kmer, _kmerSize-2-hh, binrev[NT2int(_readseq[delete_pos])]);
								}
#ifdef DEBION
								std::cout << "try corrected kmer    " << corrected_kmer.toString(_kmerSize) <<" del size" << del_size<<   std::endl;
#endif
								//check it exists, if yes correct read
								if (_bloom->contains(min(revcomp(corrected_kmer, _kmerSize), corrected_kmer)))
								{
#ifdef DEBION

									printf("apply coorec dele F  len %i pos %i\n",del_size,delete_pos);
									printf("%i %i %i\n",+delete_pos+del_size,delete_pos,_readSize-delete_pos-del_size);
									printf("old read\n%.*s\n",_readSize,_readseq);
#endif
									memmove(_readseq+delete_pos+del_size,_readseq+delete_pos,_readSize-delete_pos-del_size);
									for(int hh=0; hh<del_size; hh++)
									{
										_readseq[delete_pos+1+hh]= _readseq[delete_pos] ;
#ifdef DEBION
										printf("setting nt %i %c pos %i\n",NT2int(_readseq[delete_pos]),(_readseq[delete_pos]),delete_pos);
#endif
									}
#ifdef DEBION

									printf("new read %.*s\n",_readSize,_readseq);
#endif
									_continue_correction = true;
									_sequence->getData().setSize(_readSize);
									_local_nb_del_corrected += del_size;
									break;

								}
							}

							
						}
						

						
						
					}
					if (_continue_correction) break;

					
				}
				
				update_nb_errors_corrected(nb_errors_cor);
				
			} // end of kmers iteration over the read
			
			
		}
		
		if(nb_kmer_trusted == _readSize-_kmerSize+1){
			break;
		}
	}
	
	
//            if(! _newSeq)
//            {
//                _newSeq =  &_readseq ; //new Sequence (_readseq);
//            }
		
//            if(_use_newSeq)
//            {
//                _newSeq->setComment(_sequence->getComment());
//                _newSeq->setIndex(_sequence->getIndex());
//            }
   // printf("fin\n%i %.*s   %i  \n",(int)_newSeq->getDataSize(),(int)_newSeq->getDataSize(),_newSeq->getDataBuffer(),_newSeq->getIndex());

	
	
	Sequence * outseq;
//            if(_use_newSeq)
//                outseq = _newSeq;
//            else
		outseq =  _sequence ;
	
	#ifdef SERIAL
			_outbank->insert (*outseq); //output corrected sequence
	#else
			_bankwriter->insert (*outseq);//
			_temp_nb_seq_done ++;
			
			if(_temp_nb_seq_done == 10000) // careful with this val, should be a divisor of buffer size in ordered bank
			{
				_bankwriter->incDone(_temp_nb_seq_done);
				_temp_nb_seq_done=0;
			}
	#endif
	
}

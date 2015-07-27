#include <stdio.h>
#include <stdlib.h>
#include  <string.h>
#include <iostream>
#include <map>
#include <string>


//inf oper genome
typedef struct {
	int64_t ref_nb_bases ;
	float ref_cov;
	int64_t nb_bases_total ;
	int64_t nb_dist_total ;
} info_t;


//info per contig from tsv file
typedef struct {
	int64_t nb_bases ;
	float ref_cov ;
	std::string genome;
} info_contig;

int main(int argc, char** argv)
{
    
    
    if (argc!=3){
        printf(" Utilisation :\n");
        printf(" gener_alea <sam_file> <tsv file>\n");
        return EXIT_FAILURE;
    }
	
	
	
	std::map< std::string , info_t> infomap;
	
	info_t curinfo;
	
	FILE * sam_file =  fopen(argv[1],"r");

	
	FILE * contig_file =  fopen(argv[2],"r");
	std::map< std::string , info_contig> cmap;

	
	
	char * tempread = (char *) calloc(1000000,1);
	
	char * tempref = (char *)  malloc(1000000);
	char * ref_name = (char *)  malloc(1000000);
	
	int64_t nb_bases_total = 0;
	int64_t nb_dist_total = 0;
	int64_t read_parsed =0;
	int NM = 0;
	char *line = NULL;
	size_t linecap = 0;
	ssize_t linelen;
	
	
	
	char * tempgeno = (char *)  malloc(1000000);
	char * tempcontig = (char *)  malloc(1000000);
	float c_cov;
	int clen;
	info_contig cur_ci;
	while ((linelen = getline(&line, &linecap, contig_file)) > 0)
	{
		
		sscanf(line, "%s %s %*s %f %i %*s", tempgeno,tempcontig,&c_cov,&clen);
	
		cur_ci= cmap[std::string(tempcontig)];
		cur_ci.genome = std::string(tempgeno);
		cur_ci.ref_cov = c_cov;
		cur_ci.nb_bases = clen;
		cmap[std::string(tempcontig)] = cur_ci;
	}
	
	
	//set geno size
	for (std::map< std::string , info_contig>::iterator it=cmap.begin(); it!=cmap.end(); ++it)
	{
		cur_ci = it->second;
		
		curinfo = infomap[cur_ci.genome] ;
		curinfo.ref_nb_bases += cur_ci.nb_bases;
		infomap[cur_ci.genome]  = curinfo;
	}
	
	
	
	
	printf("Parsing sam file \n");
	while ((linelen = getline(&line, &linecap, sam_file)) > 0)
	{
		if(line[0]=='@') continue ;
		sscanf(line, "%*s %*s %s %*s %*s %*s %*s %*s %*s %s %*s %*s %*s %*s %*s %*s NM:i:%i", tempref,tempread,&NM);
		
		sscanf(tempref,"%*[^|]|%*[^|]|%*[^|]|%[^|]",ref_name);
		

		std::string refname = std::string(ref_name);
		
		std::string genoname = cmap[refname].genome;
		
		//printf("%s mapped to %s \n",ref_name,genoname.c_str());
		curinfo = infomap[genoname] ;
		curinfo.nb_bases_total +=  strlen (tempread);
		curinfo.nb_dist_total += NM ;
		curinfo.ref_cov = cmap[refname].ref_cov;
		infomap[genoname] = curinfo ;
		

		
		nb_bases_total+= strlen (tempread);
		nb_dist_total+= NM;
		read_parsed++;
		//printf("read %s  contig %s  ref %s %i \n",tempread,tempref,ref_name,NM);
		//printf(" %s \n",ref_name);
		
		if(read_parsed % 500000 ==0) {printf("%lli\t",read_parsed);fflush(stdout);}
	}
	//printf("\n");
	//	fwrite(line, linelen, 1, stdout);

	printf("NB genomes %lu \n",infomap.size());
	printf("genome\terr_rate\tnb_bases_aligned\tnb_err\tref_cov\tmeas_cov\tgeno_size\n");

	// show content:
	for (std::map< std::string,info_t>::iterator it=infomap.begin(); it!=infomap.end(); ++it)
	{
		printf("%s\t%6.3g\t%10lli\t%10lli\t%10.3f\t%10.3f\t%lli\n", it->first.c_str(),100.0* it->second.nb_dist_total/ (float) it->second.nb_bases_total,it->second.nb_bases_total,it->second.nb_dist_total,it->second.ref_cov,
			   (it->second.nb_bases_total)/((float)(it->second.ref_nb_bases)),
			   it->second.ref_nb_bases
			   );

	}
	
	
	printf("nb total bases : %lli\n",nb_bases_total);
	printf("nb total mis  : %lli   %2g\n",nb_dist_total,100.0* nb_dist_total/ (float) nb_bases_total);
	
	free(tempread);
	
	
	fclose(contig_file);
	fclose(sam_file);
}


																																								  
								

//gpl thorfinn@binf.ku.dk
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/thread_pool.h>
#include <htslib/kstring.h>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <map>
#include <vector>
#include <cmath>
#include <getopt.h>
#include <ctime>

typedef struct{
  bam1_t **d;//data
  unsigned l;//at pos
  unsigned m;//maxpos;
}queue_t;

size_t nproc=0;
int nreads_per_pos=4;//assuming this is the number reads per pos. If larger the program will reallocate
char out_mode[5]="wb";

queue_t *init_queue_t(int l){
  //  fprintf(stderr,"initializing queue with: %d elements\n",l);
  queue_t *ret =(queue_t *) malloc(sizeof(queue_t));
  ret->d =(bam1_t **) malloc(l*sizeof(bam1_t*));
  ret->l=0;
  ret->m=l;
  for(int i=0;i<ret->m;i++){
    //    fprintf(stderr,"queue[%d] init\n",i);
    ret->d[i] = bam_init1();
  }
  return ret;
}

void realloc_queue(queue_t *q){
  //  fprintf(stderr,"reallcing from: q:%d to q:%d\n",q->m,2*q->m);

  for(int i=0;0&&i<q->l;i++)
    fprintf(stderr,"inqueu[%d].pos:%d\n",i,q->d[i]->core.pos);

  bam1_t **d2 = (bam1_t **) malloc(2*q->m*sizeof(bam1_t*));

  for(int i=0;i<q->l;i++)
    d2[i] = q->d[i];
  for(int i=q->l;i<2*q->m;i++){
    d2[i] = bam_init1();
    d2[i]->core.pos=-1;
  }

  free(q->d);
  q->d=d2;
  q->m=2*q->m;
  
  for(int i=0;0&&i<q->m;i++)
    fprintf(stderr,"onqueu[%d].pos:%d\n",i,q->d[i]->core.pos);
  
}

htsFormat *dingding2 =(htsFormat*) calloc(1,sizeof(htsFormat));


char *mystr =new char[2048];

kstring_t *cigs =(kstring_t*) calloc(0,sizeof(kstring_t));
kstring_t *refs =(kstring_t*) calloc(0,sizeof(kstring_t));


int do_magic(queue_t *q,bam_hdr_t *hdr,samFile *fp){
  //fprintf(stderr,"do_magic queue->l:%d queue->m:%d chr:%d pos:%d\n",q->l,q->m,q->d[0]->core.tid,q->d[0]->core.pos);
 
  for(int i=0;i<q->l;i++)
    assert(sam_write1(fp, hdr,q->d[i])>=0);
  
  return q->l;
}

//par names are the same as global ones. Shouldnt matter, but maybe better to avoid global variables
samFile *get_outfile(char *fn_out,char out_mode[5],htsFormat *dingding2,int nthreads,htsThreadPool &p){
  samFile *out = NULL;
  static int counter = 1;//instead of global
  char fn_out2[strlen(fn_out)+10];
  snprintf(fn_out2,strlen(fn_out)+100,"%d.%s",counter,fn_out);
  if ((out = sam_open_format(fn_out2, out_mode, dingding2)) == 0) {
    fprintf(stderr,"Error opening file for writing\n");
    exit(1);
  }
  fprintf(stderr,"\t -> Has opened file: %s\n",fn_out2);

  //should this be done once or for every file
  if(counter==1&&nthreads>1){
    if (!(p.pool = hts_tpool_init(nthreads))) {
      fprintf(stderr, "Error creating thread pool\n");
      exit(2);
    }
    if (out) hts_set_opt(out, HTS_OPT_THREAD_POOL, &p);
    
  }
  counter++;
  if(counter>10)
    exit(0);
  return out;
}


int usage(FILE *fp, int is_long_help)
{
    fprintf(fp,
"\n"
"Usage: ./qname_splitter [options] <in.bam>|<in.sam>|<in.cram> \n"
"\n"
"Options:\n"
// output options
"  -b       output BAM\n"
"  -C       output CRAM (requires -T)\n"
"  -o FILE  output file name \n"
"  -T FILE  reference in the fastaformat (required from reading and writing crams)\n"
"  -@ INT   Number of threads to use\n"
	    "  -n INT   Number alignments in each output file\n"
// read filters
	    );
    fprintf(fp,
	    "\nNotes:\n");
    return 0;
}


int main(int argc, char **argv){
  clock_t t=clock();
  time_t t2=time(NULL);

  char *fname,*refName;
  samFile *in=NULL;
  samFile *out=NULL;
  fname=refName=NULL;
  char *fn_out = NULL;
  int c;
  int nthreads = 1;
  htsThreadPool p = {NULL, 0};
  int naln = 100;
  if(argc==1){
    usage(stdout,0);
    return 0;
  }
  
  static struct option lopts[] = {
    {"add", 1, 0, 0},
    {"append", 0, 0, 0},
    {"delete", 1, 0, 0},
    {"verbose", 0, 0, 0},
    {"create", 1, 0, 'c'},
    {"file", 1, 0, 0},
    {NULL, 0, NULL, 0}
  };
  
  while ((c = getopt_long(argc, argv,
			  "bCo:T:p:@:n:",
			  lopts, NULL)) >= 0) {
    switch (c) {
        case 'b': out_mode[1] = 'b'; break;
        case 'C': out_mode[1] = 'c'; break;
        case 'T': refName = strdup(optarg); break;
        case 'o': fn_out = strdup(optarg); break;
        case '@': nthreads = atoi(optarg); break;
    case 'n': naln = atoi(optarg); break;
        case '?':
	  if (optopt == '?') {  // '-?' appeared on command line
	    return usage(stdout,0);
	  } else {
	    if (optopt) { // Bad short option
	      fprintf(stdout,"./qname_splitter invalid option -- '%c'\n", optopt);
	    } else { // Bad long option
	      // Do our best.  There is no good solution to finding
	      // out what the bad option was.
	      // See, e.g. https://stackoverflow.com/questions/2723888/where-does-getopt-long-store-an-unrecognized-option
	      if (optind > 0 && strncmp(argv[optind - 1], "--", 2) == 0) {
		fprintf(stdout,"./qname_splitter unrecognised option '%s'\n",argv[optind - 1]);
	      }
	    }
	    return 0;//usage(stderr, 0);
	  }
    default:
      fprintf(stderr,"adsadsfasdf\n");
      fname = strdup(optarg);
      fprintf(stderr,"assinging: %s to fname:%s\n",optarg,fname);
      break;
    }
  }
  fname = strdup(argv[optind]);

  if(!fname){
    fprintf(stderr,"\t-> No input file specified\n");
    usage(stdout,0);
    return 0;
  }

  if(!fn_out){
    fprintf(stderr,"\t-> No output file specified\n");
    usage(stdout,0);
    return 0;
  }
  
  fprintf(stderr,"./qname_splitter refName:%s fname:%s out_mode:%s nthread:%d\n",refName,fname,out_mode,nthreads);
  
  if(refName){
    char *ref =(char*) malloc(10 + strlen(refName) + 1);
    snprintf(ref,10 + strlen(refName), "reference=%s", refName);
    hts_opt_add((hts_opt **)&dingding2->specific,ref);
    free(ref);
  }

  if(strstr(fname,".cram")!=NULL &&out_mode[1]=='c'&&refName==NULL){
    fprintf(stderr,"\t-> cram file requires reference with -T FILE.fa \n");
    return 0;
  }
  if(out_mode[1]=='c'&&refName==NULL){
    fprintf(stderr,"\t-> cram file requires reference with -T FILE.fa \n");
    return 0;
  }
  
  if((in=sam_open_format(fname,"r",dingding2))==NULL ){
    fprintf(stderr,"[%s] nonexistant file: %s\n",__FUNCTION__,fname);
    exit(0);
  }

  bam_hdr_t  *hdr = sam_hdr_read(in);  
  out = get_outfile(fn_out,out_mode,dingding2,nthreads,p);
  assert(sam_hdr_write(out, hdr) == 0);
  
  queue_t *queue = init_queue_t(nreads_per_pos);  
  

  
  bam1_t *b = bam_init1();

  int ret;
  
  // char *tmpnam =(char*) calloc(2048,sizeof(char));
  size_t nwrite = 0;
  while(((ret=sam_read1(in,hdr,b)))>0){
    nproc++;
    //catch case where there is one read in queue, and the next read is a new read
    //then we simply write the element in the queue to the output
    if(queue->l==1 && strcasecmp(bam_get_qname(queue->d[0]),bam_get_qname(b))!=0){
      nwrite++;
      assert(sam_write1(out, hdr, queue->d[0])>=0);
      queue->l =0;
    }
    if(queue->l>1 && strcasecmp(bam_get_qname(queue->d[0]),bam_get_qname(b))!=0){
      //      fprintf(stderr,"calling do_magic\n");
      nwrite += do_magic(queue,hdr,out);
      queue->l =0;
      if(nwrite>=naln){
	sam_close(out);
	out = get_outfile(fn_out,out_mode,dingding2,nthreads,p);
	assert(sam_hdr_write(out, hdr) == 0);
	nwrite = 0;
      }
    }

    if(queue->l==queue->m)
      realloc_queue(queue);

    bam_copy1(queue->d[queue->l++],b);
  }
  do_magic(queue,hdr,out);
  queue->l;
  assert(sam_close(out)==0);
  assert(sam_close(in)==0);
  for(int i=0;i<queue->m;i++)
    bam_destroy1(queue->d[i]);
  free(queue->d);
  free(queue);
  bam_hdr_destroy(hdr);

  delete [] mystr;
  bam_destroy1(b);
  hts_opt_free((hts_opt *)dingding2->specific);
  free(dingding2);
  fprintf(stderr,"    Dumpingfiles:\t\'%s\'\n",fn_out);
  free(fn_out);
  free(fname);
 
  fprintf(stderr,
	  "\t[ALL done] cpu-time used =  %.2f sec\n"
	  "\t[ALL done] walltime used =  %.2f sec\n"
	  ,(float)(clock() - t) / CLOCKS_PER_SEC, (float)(time(NULL) - t2));  

  return 0;
}

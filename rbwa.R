

PRIVATE_BWA_LIBRARIES_LOADED=FALSE;

private_load_bwa_libraries<-function()
	{
	if(!PRIVATE_BWA_LIBRARIES_LOADED)
		{
		##load.dynamic.libraries(c("libbwa.so","librbwa.so"))
		dyn.load("libbwa.so");
		dyn.load("librbwa.so");		
		}
	PRIVATE_BWA_LIBRARIES_LOADED=TRUE;
	}

#' Opens a BWA index
#'
#' @param filename the reference file
#' @keywords bwa index
#' @return bwa context
#' @examples
#' bwa.open('file.h5')
bwa.open<-function(filename)
	{
	private_load_bwa_libraries();
	.Call("RBwaOpen",filename);
	}


#' Close an ibd context and realease the associated resources
#'
#' @param bwa the BWA context
#' @keywords bwa
#' @examples
#' bwa.close(ibd)
bwa.close<-function(bwa)
	{
	private_load_bwa_libraries();
	.Call("RBwaClose",bwa)
	}

#' map one read
#'
#' @param bwa the BWA context
#' @param readseqR the read sequence
#' @keywords bwa map
#' @examples
#' ibd.close(ibd)
bwa.map<-function(bwa,readseqR)
	{
	private_load_bwa_libraries();
	.Call("RBwaMap",bwa, readseqR)
	}

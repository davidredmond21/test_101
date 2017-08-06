#include<R.h>
#include<vector>
extern "C" {
	void which_ix( double * vec, int * maxidx, int* len ) {
	// vec is the array of numbers  as input
	// location is the integer for index of max element 
	// n is the length of array
	// Operating OS is Ubuntu
	// David Redmond : id 14207377
	//
  int ilength = *len; // de reference the pointer to len
 	bool many = false;  // boolean flag to detect oif there are multiple max
	*maxidx = 1;// initialize the de-referenced "location" value 
	double temp_max = vec[0]; 	// initialize the max to the 1st element
	for ( int j = 1; j < ilength; j ++ ) { 
	  // Rprintf("detected number of elements in array %d \r\n", ilength ); 
	  if ( vec[j] > temp_max)
		{ 	temp_max = vec[j] ; 
			  *maxidx  = j+1;
		}
		else if ( vec[j]== temp_max) 
			many = true ;
		} // end for loop 
		if( many == true)
	 	Rprintf("detected multiple maxima the array the first is %d \r\n", *maxidx ); 
	} // end sumCPP
}


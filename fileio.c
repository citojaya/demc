#include "common.h"

//*--- Find the reading position in the file--*/
void FindRec(FILE *inFile, char* strDest)
{
	int nbytes = 256;
	char* strSrc;
	strSrc = (char *)malloc(nbytes+1);

	rewind(inFile);
	int n=strlen(strDest);
	while(!feof(inFile))
	{
		fgets(strSrc, 256, inFile);
		strSrc[n]='\0';
		if (strcmp(strDest, strSrc) == 0)
		{
			break;
		}
	}

	if(strcmp(strDest, strSrc) != 0)
	{
		//free(strSrc);
		//printf("Unable to find relevant info of: %s \n", strDest);
		exit(1);
	}
	free(strSrc);
}

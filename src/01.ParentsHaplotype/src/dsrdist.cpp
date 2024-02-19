#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <sstream>

int main(int argh, char** argv){
	char* ifile;
	char* ofile;
	if(argh == 3){
		ifile = argv[1];
		ofile = argv[2];
	}else{
		ifile = (char*)"3.out";
		ofile = (char*)"4.out";
	}

	int maxcliquecount = 0;
	char buffer[51200];
	buffer[51200-1] = '\0';
	std::ifstream iHMP(ifile);
	if (! iHMP.is_open()) {
		std::cerr << "Error: Invalid file " << ifile << std::endl;
		return -1;
	}

	// header
	iHMP.getline(buffer, 51200);
	int MCcount = 0;
	for(int i=0; i<51200; i++){
		if(buffer[i] == '\0'){
			if(i != 0){
				MCcount += 1;
			}
			break;
		}
		if(buffer[i] == '\t'){
			MCcount += 1;
		}
	}

	if(MCcount == 0){
		std::ofstream OFILE(ofile);
		if (! OFILE.is_open()) {
			std::cerr << "Error: Invalid file " << ifile << std::endl;
			return -1;
		}
		OFILE << "";
		OFILE.close();
		return 0;
	}

	// SNP matrix
	unsigned N_sites = 0;
	unsigned binsize = 1000000;
	unsigned** N_diff = new unsigned*[MCcount];
	for(int i=0; i<MCcount; i++){
		N_diff[i] = new unsigned[MCcount];
		for(int j=0; j<MCcount; j++){
			N_diff[i][j] = 0;
		}
	}
	unsigned** N_missing = new unsigned*[MCcount];
	for(int i=0; i<MCcount; i++){
		N_missing[i] = new unsigned[MCcount];
		for(int j=0; j<MCcount; j++){
			N_missing[i][j] = 0;
		}
	}
	unsigned char* GTlst = new unsigned char[MCcount];
	char GTbuffer[3];
	unsigned char GTbuffercursor = 0;
	int MCcursor = 0;
	while(iHMP.getline(buffer, 51200)){
		N_sites += 1;
		// init GT vector
		MCcursor = 0;
		GTbuffercursor = 0;
		for(int i=0; i<51200; i++){
			if(buffer[i] == '\0'){
				if(i>0 && buffer[i-1] != '\t'){
					if(GTbuffer[0] == '-'){
						GTlst[MCcursor] = 3;
						break;
					}else if(GTbuffer[0] == '0'){
						GTlst[MCcursor] = 0;
						break;
					}else if(GTbuffer[0] == '1'){
						GTlst[MCcursor] = 1;
						break;
					}else if(GTbuffer[0] == '2'){
						GTlst[MCcursor] = 2;
						break;
					}
				}else{
					break;
				}
			}else if(buffer[i] == '\t'){
				if(GTbuffer[0] == '-'){
					GTlst[MCcursor] = 3;
				}else if(GTbuffer[0] == '0'){
					GTlst[MCcursor] = 0;
				}else if(GTbuffer[0] == '1'){
					GTlst[MCcursor] = 1;
				}else if(GTbuffer[0] == '2'){
					GTlst[MCcursor] = 2;
				}
				GTbuffercursor = 0;
				MCcursor += 1;
			}else{
				GTbuffer[GTbuffercursor] = buffer[i];
				GTbuffercursor += 1;
			}
		}

		// For checking weather the loaded input is the same with provided input
		// for(int i=0; i<MCcount; i++){
		// 	std::cout << (int)GTlst[i] << ((i==MCcount-1)?'\n':'\t');
		// }

		for(int i=0; i<MCcount-1; i++){
			for(int j=i+1; j<MCcount; j++){
				if(GTlst[i] == 3 || GTlst[j] == 3){
					N_missing[i][j] += 1;
					N_missing[j][i] += 1;
				}else{
					if(GTlst[i] != 1 && GTlst[j] != 1 && GTlst[i] != GTlst[j]){
						N_diff[i][j] += 1;
						N_diff[j][i] += 1;
					}
				}
			}
		}
	}
	iHMP.close();

	// output
	std::ofstream OFILE(ofile);
	if (! OFILE.is_open()) {
		std::cerr << "Error: Invalid file " << ifile << std::endl;
		return -1;
	}
	for(int i=0; i<MCcount; i++){
		for(int j=0; j<MCcount; j++){
			std::stringstream obuffer;
			obuffer << std::fixed << std::setprecision(3) << (double)(N_diff[i][j]) * binsize / (double)(binsize - N_missing[i][j]);
			// Note: the original "(double)(N_diff[i][j] * binsize)" may causing overflow of int data type, because "(N_diff[i][j] * binsize)" is still int and easy to reach the maxium value of 2147483647（10^9）
			OFILE << obuffer.str() << ((j==MCcount-1)?'\n':'\t');
		}
	}
	OFILE.close();


	return 0;
}
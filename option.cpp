#include "option.h"

OptionParser::OptionParser(int argc, char ** argv){
	for (int i = 1; i < argc; i += 2){
		if (argv[i][0] == '-'){
			if (i + 1 <= argc){
				std::string tmp(argv[i]);
				paras.insert(std::make_pair(tmp, argv[i + 1]));
			}else{
				isValid = false;
				return;
			}
		} else {
			isValid = false;
			return;
		}
	}
	isValid = true;
}

std::string OptionParser::getPara(const char * str){
	std::string tmp(str);
	std::map<std::string, char *>::iterator it = paras.find(tmp);
	if (it == paras.end())
		return "";
	else {
		std::string s(it->second);
		return s;
	}
}

bool OptionParser::validCheck(){
	return isValid;
}

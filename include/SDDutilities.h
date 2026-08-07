#ifndef SDD_UTILITIES_H
#define SDD_UTILITIES_H

#include <string>
#include <vector>
#include <sstream>

inline std::string trim(const std::string& s)
{
    size_t first = s.find_first_not_of(" \t\r\n");

    if(first == std::string::npos)
        return "";

    size_t last = s.find_last_not_of(" \t\r\n");

    return s.substr(first,last-first+1);
}

inline std::vector<std::string> split(const std::string& line,char delimiter)
{
    std::vector<std::string> tokens;

    std::stringstream ss(line);

    std::string item;

    while(std::getline(ss,item,delimiter))
        tokens.push_back(trim(item));

    return tokens;
}


inline std::vector<double> parseDoubleList(const std::vector<std::string>& values)
{
    std::vector<double> output;


    for(auto& v : values)
    {
        if (v.empty())
	    continue;

	output.push_back(std::stod(v));
    }


    return output;
}

inline std::vector<int> parseIntList(const std::vector<std::string>& values)
{
    std::vector<int> output;


    for(auto& v : values)
    {
	if (v.empty())
	    continue;

        output.push_back(std::stoi(v));
    }


    return output;
}



#endif

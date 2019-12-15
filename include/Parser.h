#pragma once

#include <fstream>
#include <string>

/*!
 * \class Parser
 * \brief Pass a config file to it so it parses and set input values
 */
class Parser
{
private:

public:
	/*!
	 * \fn Parser(const char* configFile)
	 * \brief Constructor
	 * \param configFile name of the config file
	 */
	Parser(const char* configFile);
	~Parser(){}

	/*!
	 * \fn void readLines(std::fstream& fs)
	 * \brief Read the config file line by line
	 * \param fs the file stream object for the config file
	 */
	void readLines(std::ifstream& fs);
};

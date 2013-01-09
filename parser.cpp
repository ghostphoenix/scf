// parser.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iterator>
#include <iostream>
#include <fstream>
using std::string;
#include <sstream>
#include <vector>
#include <string>
#include <list>
#include <algorithm>
#include <stdlib.h>
#include <map>

namespace CONSTANT
{
    const std::string SYMBOL[37] = 
    {
        "e",
        "H",                                                                                                 "He", 
        "Li", "Be",                                                            "B" , "C",  "N",  "O",  "F",  "Ne",
        "Na", "Mg",                                                            "Al", "Si", "P ", "S ", "Cl", "Ar",
        "K",  "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr"

    };    

    std::map<const std::string,int> SYMBOLMAP;

    void buildSymbolMap()
    {
        for (int i = 0; i <= 19; i++)
        {   
            SYMBOLMAP[SYMBOL[i]] = i;
        }
    }
} 

class geometry
{
private:
    int               noAtoms;
    std::vector<double>    coordinates;
    std::vector<int>       Z;

public:
    geometry(std::vector<int> &_z, std::vector<double> &c)
    {
        noAtoms     = c.size();
        Z           = _z;
        coordinates = c;
    }
};


namespace stringToolkit
{
    void trimSpaces(std::string &s)
    {
        s.erase(std::remove_if(s.begin(), s.end(), isspace), s.end());
    }
    
    void toUpper(std::string &s)
    {
        std::transform(s.begin(), s.end(),s.begin(), ::toupper);
    }

    void toUpperLower(std::string &s)
    {
        std::transform(s.begin(), s.begin(),s.begin(), ::toupper);
        if (s.size() > 1)
        std::transform(s.begin()+1, s.end(),s.begin()+1, ::tolower);
    }

    std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) 
    {
        std::stringstream ss(s);
        std::string item;
        while(std::getline(ss, item, delim)) 
        {
            elems.push_back(item);
        }
        return elems;
    }

    std::vector<std::string> split(const std::string &s, char delim) 
    {
        std::vector<std::string> elems;
        return split(s, delim, elems);
    }
}


class CSVRow
{    
    private:
        std::vector<std::string>    m_data;                

    public:
        std::string const& operator[](std::size_t index) const
        {
            return m_data[index];
        }
        std::size_t size() const
        {
            return m_data.size();
        }
        void readNextRow(std::istream& str)
        {
            std::string         line;
            std::getline(str,line);

            std::stringstream   lineStream(line);
            std::string         cell;

            m_data.clear();
            while(std::getline(lineStream,cell,','))
            {
                m_data.push_back(cell);
            }
        }

        std::string get(const int i)
        {
            return m_data[i];
        }

        double getDouble(const int i)
        {
            return atof( m_data[i].c_str() );
        }

       std::vector<double> getVectorNoZero(const int start, int end)
       {
            std::vector<double> ret;

            if (end==0) end = m_data.size() - 1;

            for (int i = start; i <= end; i++) 
            {
                double temp = this->getDouble(i);
                
                if (temp != 0.0)
                   {
                    ret.push_back(temp);
                   }
            }  
            return ret;
       }     

       int size()
       {
            return m_data.size();
       }     

       int size(int i)
       {
            return m_data.size();
       }

};

std::istream& operator>>(std::istream& str, CSVRow& data)
{
    data.readNextRow(str);
    return str;
}  

class CSVTable
{
    private:
        std::vector<CSVRow>    m_data; 
        CSVRow                 workingRow; 
        string                 filename;
        std::vector<double>    workingMatrix;
        std::vector<int>       workingMatrixInt;

    public:
        CSVRow &getRow(const int i)                     { return m_data[i]; }
        std::string getCell(const int i, const int j)   { return m_data[i].get(j); }
        void   setRow(CSVRow row)                       { m_data.push_back(row); }
        int    getNumberRows()                          { return m_data.size(); }
        
        int id;

        CSVTable(std::string fname) 
        {
            std::ifstream       file(fname);
           
            while(file >> workingRow)
            {
                this->setRow(workingRow);
            }
        }


        std::vector<double> &getMatrixDouble()
        {
            return workingMatrix;
        }

        std::vector<double> matrix()
        {
            workingMatrix.clear();
            int N = getNumberRows();
            int M;
            
            for (int i = 0; i < N; i++)
            {
                workingRow = this->getRow(i);
                M = workingRow.size(i);
                for (int j = 0; j < M; j++)
                {
                    workingMatrix.push_back(workingRow.getDouble(j));
                }
            }
            return workingMatrix;
        }


        geometry parseGeo()
        {
            workingMatrix.clear();
            workingMatrixInt.clear();
            std::cout << "\n No of Atoms :" << std::endl;
            int N = getNumberRows();
            unsigned int i, j;
            std::cout << "\n No of Atoms :" << N  << std::endl;

            for (i = 0; i < m_data.size(); i++)
            {
                workingRow = this->getRow(i);
                
                std::cout << "\nAtom 1 is :" << workingRow.get(0)  << std::endl;
                workingMatrixInt.push_back(-99);
                for (j = 1; j <= 3; j++)
                {
                    workingMatrix.push_back(workingRow.getDouble(j));
                    std::cout << workingRow.getDouble(j) <<  ", " << std::endl;
                }
            }
            
            geometry geo(workingMatrixInt,workingMatrix);
            return geo;
        }

};

void parseBasis(CSVTable &table)
{
    CONSTANT::buildSymbolMap();
    unsigned int i;
    int Z;
    int basisNo = -1;

    std::string symbol;
    std::string tmpStr;
    std::vector< std::vector<double> > exp;
    std::vector< std::vector<double> > coeff;
    std::vector<double>           tempVec;
    std::vector<std::string> temp;

 for (int iRow = 0; iRow < table.getNumberRows(); iRow++)
    {
        char first = table.getCell(iRow,0)[0];
            
        switch (first) 
        {
            case 's':
            case 'p':
            case 'd':
            case 'f':
            case 'g':
            case 'h':
                basisNo++;
                tempVec.clear();
 
                symbol = table.getCell(iRow,1);
                stringToolkit::trimSpaces(symbol);
                stringToolkit::toUpperLower(symbol);
                Z = CONSTANT::SYMBOLMAP[symbol];
                

                std::cout << "Angular momentum is :" << first  << std::endl;
                std::cout << "Element             :" << symbol << std::endl;
                std::cout << "Atomic Number =     :" << Z << std::endl;

                tempVec = table.getRow(iRow).getVectorNoZero(2,0);

                exp.push_back(tempVec);

                for (i=0; i < exp[basisNo].size(); i++) 
                {
                    std::cout << exp[basisNo][i] << ", "; 
                }    
                    std::cout << "\n";      
                break;
            case 'c':
//                tmpStr = row[2];
//                tmpStr.erase( remove( tmpStr.begin(), tmpStr.end(), ' ' ), tmpStr.end() );
//                testCon = row[2].c_str();

                temp = stringToolkit::split(table.getCell(iRow,1),'.');

                std::cout << "Contraction: " << temp[0] << '-' << temp[1] << std::endl;

                tempVec = table.getRow(iRow).getVectorNoZero(2,0);

   
                coeff.push_back(tempVec);

                for (i=0; i < coeff[basisNo].size(); i++) 
                {
                    std::cout << coeff[basisNo][i] << ", "; 
                }    
                    std::cout << "\n"; 
                break;
            case '!':
                std::cout <<  table.getCell(iRow,0) << std::endl;
                break;
            default:
                std::cout <<  table.getCell(iRow,0) << " is not recognised \n";
                break;
            }
    }
}

int _tmain(int argc, _TCHAR* argv[])
{

    CSVTable table("testBasis.csv");
    CSVTable tableMat("testMat.csv");
    CSVTable tableGeo("testGeo.csv");

    parseBasis(table);

    std::vector<double>           tempVec = tableMat.matrix();
    
    std::cout << "/n";
    for (unsigned int i = 0; i < tempVec.size(); i++)
    {
        std::cout <<  tempVec[i] << ", ";
    }

    tableGeo.parseGeo();
    return 0;
}

    


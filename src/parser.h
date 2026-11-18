/*
Expression parser
lang:   C++ code, project developed with wxDev-C++ (rewritten for ANSI C++ under MS VC 2008 and gcc 4.3 by Jörg Schulla)
author: Jos de Jong, 2007-12-22
		changes made by Jörg Schulla 2011-03-05:
		- putting all files together into one source and one header file (for easier inclusion into WaSiM)
		- removing namespace usage in headers (is generally not recommended to use namespaces in headers)
		- replacing snprintf functions by sprintf
		- replacing vsnprintf by sprintf
		- replacing error constructor with variable arg list by seperate calls with a clear argument list
		- creating messages by using stringstream instead of sprintf
		- made variablelist an external vector which can be initialized before and used after the parser goes out of life
		- removed row() from parser (we use single line expressions only in the control file, even if they are splitted into multiple rows
site:   www.speqmath.com

Features:
    Operators:
        & | << >>
        = <> < > <= >=
        + -
        * / % ||
        ^
        !

    Functions:
        Abs, Exp, Sign, Sqrt, Log, Log10
        Sin, Cos, Tan, ASin, ACos, ATan
        Factorial

    Variables:
        PI, EULER
        you can define your own variables
        
    Other:
        Scientific notation supported
        Error handling supported

*/


#ifndef PARSER_H
#define PARSER_H

// declarations
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <cmath>
#include <vector>
#include <cstdarg>

const int NAME_LEN_MAX = 30;
const int EXPR_LEN_MAX = 255;
const int ERR_LEN_MAX = 255;

void toupper(char upper[], const char str[]);
double factorial(double value);
double sign(double value);

class Error {
    public:
        Error(const int col, const int id, const double darg);
        Error(const int col, const int id, const int iarg);
        Error(const int col, const int id, const char *carg);
        Error(const int col, const int id);

        int get_col() {return err_col;} // Returns the column of the error
        int get_id() {return err_id;}   // Returns the id of the error
        char* get_msg() {return msg;}   // Returns a pointer to the error msg

    private:
        int err_col;    // column (position) where the error occured
        int err_id;     // id of the error
        char msg[255];
        
        char* msgdesc(const int id);
};



class Variablelist {
    public:
        bool exist(const char* name);
        bool add(const char* name, double value);
        bool del(const char* name);
        void clear();

        bool get_value(const char* name, double* value);
        bool get_value(const int id, double* value);
        int  get_id(const char* name);
        bool set_value(const char* name, const double value);

    private:
        struct VAR
        {
            char name[NAME_LEN_MAX+1];
            double value;
        };

		std::vector<VAR> var;
};

class Parser
{
    // public functions
    public:
        Parser(Variablelist *varlist);
        char* parse(const char expr[]);

    // enumerations
    private:
    
        enum TOKENTYPE {NOTHING = -1, DELIMETER, NUMBER, VARIABLE, FUNCTION, UNKNOWN};
    
        enum OPERATOR_ID {AND, OR, BITSHIFTLEFT, BITSHIFTRIGHT,                 // level 2
                       EQUAL, UNEQUAL, SMALLER, LARGER, SMALLEREQ, LARGEREQ,    // level 3
                       PLUS, MINUS,                     // level 4
                       MULTIPLY, DIVIDE, MODULUS, XOR,  // level 5
                       POW,                             // level 6
                       FACTORIAL};                      // level 7

    // data
    private:
        char expr[EXPR_LEN_MAX+1];    // holds the expression
        char* e;                      // points to a character in expr
        
        char token[NAME_LEN_MAX+1];   // holds the token
        TOKENTYPE token_type;         // type of the token

        double ans;                   // holds the result of the expression
        char ans_str[255];            // holds a string containing the result 
                                      // of the expression

        Variablelist *user_var;        // list with variables defined by user

    // private functions
    private:
        void getToken();
        
        double parse_level1();
        double parse_level2();
        double parse_level3();
        double parse_level4();
        double parse_level5();
        double parse_level6();
        double parse_level7();
        double parse_level8();
        double parse_level9();
        double parse_level10();
        double parse_number();

        int get_operator_id(const char op_name[]);
        double eval_operator(const int op_id, const double &lhs, const double &rhs);
        double eval_function(const char fn_name[], const double &value);
        double eval_variable(const char var_name[]);
        int col();
};

#endif

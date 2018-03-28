/* 
 * File:   evalexpression.h
 * Author: nash
 *
 * Created on March 12, 2016, 1:16 PM
 */

#ifndef EVALEXPRESSION_H
#define	EVALEXPRESSION_H

#include <string>
#include "exprtk.hpp"

using namespace std;

template <typename T>
void evaluate(string expression_string, double xvalue, double &W)
{
   typedef exprtk::symbol_table<T> symbol_table_t;
   typedef exprtk::expression<T>     expression_t;
   typedef exprtk::parser<T>             parser_t;
   
   T x = T(xvalue);
   
   symbol_table_t symbol_table;
   symbol_table.add_variable("x",x);
   symbol_table.add_constants();
   
   expression_t expression;
   expression.register_symbol_table(symbol_table);

   parser_t parser;
   parser.compile(expression_string,expression);
   
   T result = expression.value();
   
   
   
   W=result;
}

template <typename T>
double evalreturn(string expression_string, double xvalue)
{
   typedef exprtk::symbol_table<T> symbol_table_t;
   typedef exprtk::expression<T>     expression_t;
   typedef exprtk::parser<T>             parser_t;
   
   T x = T(xvalue);
   
   symbol_table_t symbol_table;
   symbol_table.add_variable("x",x);
   symbol_table.add_constants();
   
   expression_t expression;
   expression.register_symbol_table(symbol_table);

   parser_t parser;
   parser.compile(expression_string,expression);
   
   T result = expression.value();
   
   return(result);
}

#endif	/* EVALEXPRESSION_H */


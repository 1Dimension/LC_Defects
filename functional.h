/**
 * @file   functional.h
 * @author onedimension <onedimension@onedimension-PC>
 * @date   Sun Oct 26 14:34:48 2014
 * 
 * @brief  头文件，定义计算能量泛函的函数
 * 
 * 
 */

#ifndef _FUNCTIONAL_H
#define _FUNCTIONAL_H

using namespace std;

#include<iostream>

/** 
 * @brief 根据Q(r_i, \theta_j, \phi_k)计算Bulk Energy density
 * 
 * @param qijk 
 * @param fbulk 
 */
void calc_bulkenergy_landau_de(double *qijk,double *fbulk);

void calc_grad_bulkenergy_landau_de(double *qijk,double *grad_fbulk);

double calc_energy_Fb(double *qnlm,double *qijk,double *fbulk); ///计算bulk energy








#endif

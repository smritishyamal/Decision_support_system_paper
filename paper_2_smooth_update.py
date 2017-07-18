# -*- coding: utf-8 -*-
"""
Created on Sat Jul 30 00:34:26 2016

@author: User
"""

#! ==================================================================
from __future__ import division
from casadi import *
from numpy import *
from pylab import *
import numpy as NP
import matplotlib.pyplot as plt
from scipy import linalg, matrix
import time
NP.random.seed(0)

#%%################################################
# EAF model deatils
################################################
dae = DaeBuilder()

# State variable declarations
sm_b_0 = dae.add_x('sm_b_0')
sm_b_1 = dae.add_x('sm_b_1')
sm_b_2 = dae.add_x('sm_b_2')
sm_b_3 = dae.add_x('sm_b_3')
sm_b_4 = dae.add_x('sm_b_4')
sm_b_5 = dae.add_x('sm_b_5')
sm_b_6 = dae.add_x('sm_b_6')
sm_b_7 = dae.add_x('sm_b_7')
sm_m_cao = dae.add_x('sm_m_cao')
sm_m_c_float = dae.add_x('sm_m_c_float')
sm_m_dol = dae.add_x('sm_m_dol')
sm_e = dae.add_x('sm_e')
ss_m_ss = dae.add_x('ss_m_ss')
ss_t = dae.add_x('ss_t')
mm_n_0 = dae.add_x('mm_n_0')
mm_n_1 = dae.add_x('mm_n_1')
mm_n_2 = dae.add_x('mm_n_2')
mm_n_3 = dae.add_x('mm_n_3')
mm_n_4 = dae.add_x('mm_n_4')
mm_t = dae.add_x('mm_t')
gs_b_0 = dae.add_x('gs_b_0')
gs_b_1 = dae.add_x('gs_b_1')
gs_b_2 = dae.add_x('gs_b_2')
gs_b_3 = dae.add_x('gs_b_3')
gs_n_oil_gas = dae.add_x('gs_n_oil_gas')
gs_t = dae.add_x('gs_t')
rd_t_roof = dae.add_x('rd_t_roof')
rd_t_wall = dae.add_x('rd_t_wall')

# Add the disturbance state
disturb_add = 1
if disturb_add == 1:
    d_sm_b_3 = dae.add_x('d_sm_b_3') # Disturbance state for sm_b_3
                
# Algebraic variable declarations
gs_sum_lambda_0 = dae.add_z('gs_sum_lambda_0')
gs_sum_lambda_1 = dae.add_z('gs_sum_lambda_1')
gs_sum_lambda_2 = dae.add_z('gs_sum_lambda_2')
gs_sum_lambda_3 = dae.add_z('gs_sum_lambda_3')
gs_sum_lambda_4 = dae.add_z('gs_sum_lambda_4')
gs_sum_lambda_5 = dae.add_z('gs_sum_lambda_5')
gs_sum_lambda_6 = dae.add_z('gs_sum_lambda_6')
gs_dg_0 = dae.add_z('gs_dg_0')
gs_dg_1 = dae.add_z('gs_dg_1')
gs_dg_2 = dae.add_z('gs_dg_2')
gs_dg_3 = dae.add_z('gs_dg_3')
gs_dg_4 = dae.add_z('gs_dg_4')
gs_dg_5 = dae.add_z('gs_dg_5')
gs_dg_6 = dae.add_z('gs_dg_6')
gs_xl_0 = dae.add_z('gs_xl_0')
gs_xl_1 = dae.add_z('gs_xl_1')
gs_xl_2 = dae.add_z('gs_xl_2')
gs_xl_3 = dae.add_z('gs_xl_3')
gs_xl_4 = dae.add_z('gs_xl_4')
gs_xl_5 = dae.add_z('gs_xl_5')
gs_xl_6 = dae.add_z('gs_xl_6')
gs_xl_7 = dae.add_z('gs_xl_7')
gs_lambda_c = dae.add_z('gs_lambda_c')
gs_lambda_o = dae.add_z('gs_lambda_o')
gs_lambda_h = dae.add_z('gs_lambda_h')
gs_ntx = dae.add_z('gs_ntx')
gs_yx_0 = dae.add_z('gs_yx_0')
gs_yx_1 = dae.add_z('gs_yx_1')
gs_yx_2 = dae.add_z('gs_yx_2')
gs_yx_3 = dae.add_z('gs_yx_3')
gs_yx_4 = dae.add_z('gs_yx_4')
gs_yx_5 = dae.add_z('gs_yx_5')
gs_yx_6 = dae.add_z('gs_yx_6')
gs_yx_7 = dae.add_z('gs_yx_7')
gs_y_dry_0 = dae.add_z('gs_y_dry_0')
gs_y_dry_1 = dae.add_z('gs_y_dry_1')
gs_y_dry_2 = dae.add_z('gs_y_dry_2')
gs_y_dry_3 = dae.add_z('gs_y_dry_3')
gs_y_dry_4 = dae.add_z('gs_y_dry_4')
gs_y_dry_5 = dae.add_z('gs_y_dry_5')
gs_y_dry_6 = dae.add_z('gs_y_dry_6')
gs_y_dry_7 = dae.add_z('gs_y_dry_7')
gs_f_in_nm3_2 = dae.add_z('gs_f_in_nm3_2')
gs_f_in_nm3_3 = dae.add_z('gs_f_in_nm3_3')
gs_f_in_nm3_5 = dae.add_z('gs_f_in_nm3_5')
gs_f_in_nm3_6 = dae.add_z('gs_f_in_nm3_6')
gs_n_add_0 = dae.add_z('gs_n_add_0')
gs_n_add_1 = dae.add_z('gs_n_add_1')
gs_n_add_2 = dae.add_z('gs_n_add_2')
gs_n_add_3 = dae.add_z('gs_n_add_3')
gs_n_add_4 = dae.add_z('gs_n_add_4')
gs_n_add_5 = dae.add_z('gs_n_add_5')
gs_n_add_6 = dae.add_z('gs_n_add_6')
gs_n_add_7 = dae.add_z('gs_n_add_7')
gs_n_sm_0 = dae.add_z('gs_n_sm_0')
gs_n_sm_2 = dae.add_z('gs_n_sm_2')
gs_n_ngrs_0 = dae.add_z('gs_n_ngrs_0')
gs_n_ngrs_1 = dae.add_z('gs_n_ngrs_1')
gs_n_ngrs_2 = dae.add_z('gs_n_ngrs_2')
gs_n_ngrs_3 = dae.add_z('gs_n_ngrs_3')
gs_n_ngrs_4 = dae.add_z('gs_n_ngrs_4')
gs_n_ngrs_5 = dae.add_z('gs_n_ngrs_5')
gs_n_ngrs_6 = dae.add_z('gs_n_ngrs_6')
gs_n_ngrs_7 = dae.add_z('gs_n_ngrs_7')
gs_n_in_0 = dae.add_z('gs_n_in_0')
gs_n_in_1 = dae.add_z('gs_n_in_1')
gs_n_in_2 = dae.add_z('gs_n_in_2')
gs_n_in_3 = dae.add_z('gs_n_in_3')
gs_n_in_4 = dae.add_z('gs_n_in_4')
gs_n_in_5 = dae.add_z('gs_n_in_5')
gs_n_in_6 = dae.add_z('gs_n_in_6')
gs_n_in_7 = dae.add_z('gs_n_in_7')
gs_n_ex_0 = dae.add_z('gs_n_ex_0')
gs_n_ex_1 = dae.add_z('gs_n_ex_1')
gs_n_ex_2 = dae.add_z('gs_n_ex_2')
gs_n_ex_3 = dae.add_z('gs_n_ex_3')
gs_n_ex_4 = dae.add_z('gs_n_ex_4')
gs_n_ex_5 = dae.add_z('gs_n_ex_5')
gs_n_ex_7 = dae.add_z('gs_n_ex_7')
gs_n_oil_in = dae.add_z('gs_n_oil_in')
gs_n_push = dae.add_z('gs_n_push')
gs_n_suck = dae.add_z('gs_n_suck')
gs_n_net = dae.add_z('gs_n_net')
gs_t_ss = dae.add_z('gs_t_ss')
gs_dum7 = dae.add_z('gs_dum7')
gs_dum8 = dae.add_z('gs_dum8')
gs_f_in_o2_oxy_burn = dae.add_z('gs_f_in_o2_oxy_burn')
gs_o2_switch1 = dae.add_z('gs_o2_switch1')
gs_o2_switch2 = dae.add_z('gs_o2_switch2')
gs_o2_switch3 = dae.add_z('gs_o2_switch3')
gs_f_in_o2_jetbox1 = dae.add_z('gs_f_in_o2_jetbox1')
gs_f_in_o2_jetbox2 = dae.add_z('gs_f_in_o2_jetbox2')
gs_f_in_o2_jetbox3 = dae.add_z('gs_f_in_o2_jetbox3')
gs_jetgas1 = dae.add_z('gs_jetgas1')
gs_jetgas2 = dae.add_z('gs_jetgas2')
gs_jetgas3 = dae.add_z('gs_jetgas3')
gs_tt = dae.add_z('gs_tt')
gs_ttsm = dae.add_z('gs_ttsm')
gs_cpdt_0 = dae.add_z('gs_cpdt_0')
gs_cpdt_1 = dae.add_z('gs_cpdt_1')
gs_cpdt_2 = dae.add_z('gs_cpdt_2')
gs_cpdt_3 = dae.add_z('gs_cpdt_3')
gs_cpdt_4 = dae.add_z('gs_cpdt_4')
gs_cpdt_5 = dae.add_z('gs_cpdt_5')
gs_cpdt_6 = dae.add_z('gs_cpdt_6')
gs_cpdt_7 = dae.add_z('gs_cpdt_7')
gs_cptdt_0 = dae.add_z('gs_cptdt_0')
gs_cptdt_1 = dae.add_z('gs_cptdt_1')
gs_cptdt_2 = dae.add_z('gs_cptdt_2')
gs_cptdt_3 = dae.add_z('gs_cptdt_3')
gs_cptdt_4 = dae.add_z('gs_cptdt_4')
gs_cptdt_5 = dae.add_z('gs_cptdt_5')
gs_cptdt_6 = dae.add_z('gs_cptdt_6')
gs_dumgas = dae.add_z('gs_dumgas')
gs_dumgas3 = dae.add_z('gs_dumgas3')
gs_t_sm = dae.add_z('gs_t_sm')
gs_cpdt_sm_0 = dae.add_z('gs_cpdt_sm_0')
gs_cpdt_sm_1 = dae.add_z('gs_cpdt_sm_1')
gs_cpdt_sm_2 = dae.add_z('gs_cpdt_sm_2')
gs_cpdt_sm_3 = dae.add_z('gs_cpdt_sm_3')
gs_cpdt_sm_4 = dae.add_z('gs_cpdt_sm_4')
gs_cpdt_sm_5 = dae.add_z('gs_cpdt_sm_5')
gs_cpdt_sm_6 = dae.add_z('gs_cpdt_sm_6')
gs_cpdt_sm_7 = dae.add_z('gs_cpdt_sm_7')
gs_fh_in = dae.add_z('gs_fh_in')
gs_fh_out = dae.add_z('gs_fh_out')
gs_q_ss = dae.add_z('gs_q_ss')
gs_q = dae.add_z('gs_q')
gs_q_cool = dae.add_z('gs_q_cool')
gs_h_sm_0 = dae.add_z('gs_h_sm_0')
gs_h_sm_1 = dae.add_z('gs_h_sm_1')
gs_h_sm_2 = dae.add_z('gs_h_sm_2')
gs_h_sm_3 = dae.add_z('gs_h_sm_3')
gs_h_sm_4 = dae.add_z('gs_h_sm_4')
gs_h_sm_5 = dae.add_z('gs_h_sm_5')
gs_h_sm_6 = dae.add_z('gs_h_sm_6')
gs_h_sm_7 = dae.add_z('gs_h_sm_7')
gs_h_0 = dae.add_z('gs_h_0')
gs_h_1 = dae.add_z('gs_h_1')
gs_h_2 = dae.add_z('gs_h_2')
gs_h_3 = dae.add_z('gs_h_3')
gs_h_4 = dae.add_z('gs_h_4')
gs_h_5 = dae.add_z('gs_h_5')
gs_h_6 = dae.add_z('gs_h_6')
gs_h_7 = dae.add_z('gs_h_7')
mm_melt_0 = dae.add_z('mm_melt_0')
mm_melt_1 = dae.add_z('mm_melt_1')
mm_melt_2 = dae.add_z('mm_melt_2')
mm_melt_3 = dae.add_z('mm_melt_3')
mm_melt_4 = dae.add_z('mm_melt_4')
mm_melt_5 = dae.add_z('mm_melt_5')
mm_f_in = dae.add_z('mm_f_in')
mm_f_ex_0 = dae.add_z('mm_f_ex_0')
mm_f_ex_1 = dae.add_z('mm_f_ex_1')
mm_f_ex_2 = dae.add_z('mm_f_ex_2')
mm_f_ex_3 = dae.add_z('mm_f_ex_3')
mm_f_ex_4 = dae.add_z('mm_f_ex_4')
mm_m_mm_0 = dae.add_z('mm_m_mm_0')
mm_m_mm_1 = dae.add_z('mm_m_mm_1')
mm_m_mm_2 = dae.add_z('mm_m_mm_2')
mm_m_mm_3 = dae.add_z('mm_m_mm_3')
mm_m_mm_4 = dae.add_z('mm_m_mm_4')
mm_m_mm_5 = dae.add_z('mm_m_mm_5')
mm_f_sm_0 = dae.add_z('mm_f_sm_0')
mm_f_sm_1 = dae.add_z('mm_f_sm_1')
mm_f_sm_2 = dae.add_z('mm_f_sm_2')
mm_f_sm_3 = dae.add_z('mm_f_sm_3')
mm_f_sm_4 = dae.add_z('mm_f_sm_4')
mm_f_sm_5 = dae.add_z('mm_f_sm_5')
mm_y_0 = dae.add_z('mm_y_0')
mm_y_1 = dae.add_z('mm_y_1')
mm_y_2 = dae.add_z('mm_y_2')
mm_y_3 = dae.add_z('mm_y_3')
mm_y_4 = dae.add_z('mm_y_4')
mm_y_5 = dae.add_z('mm_y_5')
mm_ntot = dae.add_z('mm_ntot')
mm_p_power = dae.add_z('mm_p_power')
mm_q_mm_sm = dae.add_z('mm_q_mm_sm')
mm_q_cool = dae.add_z('mm_q_cool')
mm_q_mm_ss = dae.add_z('mm_q_mm_ss')
mm_f_in_kgm_c = dae.add_z('mm_f_in_kgm_c')
mm_m_slag = dae.add_z('mm_m_slag')
mm_y_slag_0 = dae.add_z('mm_y_slag_0')
mm_y_slag_1 = dae.add_z('mm_y_slag_1')
mm_y_slag_2 = dae.add_z('mm_y_slag_2')
mm_y_slag_3 = dae.add_z('mm_y_slag_3')
mm_y_slag_4 = dae.add_z('mm_y_slag_4')
mm_y_slag_5 = dae.add_z('mm_y_slag_5')
mm_p_arc = dae.add_z('mm_p_arc')
mm_t_sm = dae.add_z('mm_t_sm')
mm_f_in_kgm_c_lance = dae.add_z('mm_f_in_kgm_c_lance')
mm_f_melt_ss = dae.add_z('mm_f_melt_ss')
mm_m_solid = dae.add_z('mm_m_solid')
mm_moltendum5 = dae.add_z('mm_moltendum5')
mm_mtot = dae.add_z('mm_mtot')
mm_e_f = dae.add_z('mm_e_f')
mm_k_ml = dae.add_z('mm_k_ml')
mm_f_o2_lnc = dae.add_z('mm_f_o2_lnc')
mm_q_rad = dae.add_z('mm_q_rad')
mm_q_wall = dae.add_z('mm_q_wall')
rd_hw_exposedwall = dae.add_z('rd_hw_exposedwall')
rd_a_wall = dae.add_z('rd_a_wall')
rd_v_scrap = dae.add_z('rd_v_scrap')
rd_m_solid = dae.add_z('rd_m_solid')
rd_hs_scrappile = dae.add_z('rd_hs_scrappile')
rd_v_bath = dae.add_z('rd_v_bath')
rd_h_bath = dae.add_z('rd_h_bath')
rd_q_gascool_0 = dae.add_z('rd_q_gascool_0')
rd_q_gascool_1 = dae.add_z('rd_q_gascool_1')
rd_q_cw_0 = dae.add_z('rd_q_cw_0')
rd_q_cw_1 = dae.add_z('rd_q_cw_1')
rd_ua_t2 = dae.add_z('rd_ua_t2')
rd_m_liquid = dae.add_z('rd_m_liquid')
rd_t_gs = dae.add_z('rd_t_gs')
rd_e_f = dae.add_z('rd_e_f')
rd_q_wall = dae.add_z('rd_q_wall')
rd_q_radroof = dae.add_z('rd_q_radroof')
rd_q_radwall = dae.add_z('rd_q_radwall')
sm_sum_lambda_0 = dae.add_z('sm_sum_lambda_0')
sm_sum_lambda_1 = dae.add_z('sm_sum_lambda_1')
sm_sum_lambda_2 = dae.add_z('sm_sum_lambda_2')
sm_sum_lambda_3 = dae.add_z('sm_sum_lambda_3')
sm_sum_lambda_4 = dae.add_z('sm_sum_lambda_4')
sm_sum_lambda_5 = dae.add_z('sm_sum_lambda_5')
sm_sum_lambda_6 = dae.add_z('sm_sum_lambda_6')
sm_sum_lambda_7 = dae.add_z('sm_sum_lambda_7')
sm_sum_lambda_8 = dae.add_z('sm_sum_lambda_8')
sm_sum_lambda_9 = dae.add_z('sm_sum_lambda_9')
sm_sum_lambda_10 = dae.add_z('sm_sum_lambda_10')
sm_sum_lambda_11 = dae.add_z('sm_sum_lambda_11')
sm_sum_lambda_12 = dae.add_z('sm_sum_lambda_12')
sm_sum_lambda_13 = dae.add_z('sm_sum_lambda_13')
sm_lambda_c = dae.add_z('sm_lambda_c')
sm_lambda_o = dae.add_z('sm_lambda_o')
sm_lambda_fe = dae.add_z('sm_lambda_fe')
sm_lambda_mg = dae.add_z('sm_lambda_mg')
sm_lambda_mn = dae.add_z('sm_lambda_mn')
sm_lambda_si = dae.add_z('sm_lambda_si')
sm_lambda_al = dae.add_z('sm_lambda_al')
sm_dg_0 = dae.add_z('sm_dg_0')
sm_dg_1 = dae.add_z('sm_dg_1')
sm_dg_2 = dae.add_z('sm_dg_2')
sm_dg_3 = dae.add_z('sm_dg_3')
sm_dg_4 = dae.add_z('sm_dg_4')
sm_dg_5 = dae.add_z('sm_dg_5')
sm_dg_6 = dae.add_z('sm_dg_6')
sm_dg_7 = dae.add_z('sm_dg_7')
sm_dg_8 = dae.add_z('sm_dg_8')
sm_dg_9 = dae.add_z('sm_dg_9')
sm_dg_10 = dae.add_z('sm_dg_10')
sm_dg_11 = dae.add_z('sm_dg_11')
sm_dg_12 = dae.add_z('sm_dg_12')
sm_dg_13 = dae.add_z('sm_dg_13')
sm_xl_0 = dae.add_z('sm_xl_0')
sm_xl_1 = dae.add_z('sm_xl_1')
sm_xl_2 = dae.add_z('sm_xl_2')
sm_xl_3 = dae.add_z('sm_xl_3')
sm_xl_4 = dae.add_z('sm_xl_4')
sm_xl_5 = dae.add_z('sm_xl_5')
sm_xl_6 = dae.add_z('sm_xl_6')
sm_xl_7 = dae.add_z('sm_xl_7')
sm_xl_8 = dae.add_z('sm_xl_8')
sm_xl_9 = dae.add_z('sm_xl_9')
sm_xl_10 = dae.add_z('sm_xl_10')
sm_xl_11 = dae.add_z('sm_xl_11')
sm_xl_12 = dae.add_z('sm_xl_12')
sm_xl_13 = dae.add_z('sm_xl_13')
sm_xl_14 = dae.add_z('sm_xl_14')
sm_xl_15 = dae.add_z('sm_xl_15')
sm_t = dae.add_z('sm_t')
sm_f_in_0 = dae.add_z('sm_f_in_0')
sm_f_in_1 = dae.add_z('sm_f_in_1')
sm_f_in_2 = dae.add_z('sm_f_in_2')
sm_f_in_3 = dae.add_z('sm_f_in_3')
sm_f_in_4 = dae.add_z('sm_f_in_4')
sm_f_in_5 = dae.add_z('sm_f_in_5')
sm_f_in_6 = dae.add_z('sm_f_in_6')
sm_f_in_7 = dae.add_z('sm_f_in_7')
sm_f_in_8 = dae.add_z('sm_f_in_8')
sm_f_in_9 = dae.add_z('sm_f_in_9')
sm_f_in_10 = dae.add_z('sm_f_in_10')
sm_f_in_11 = dae.add_z('sm_f_in_11')
sm_f_in_12 = dae.add_z('sm_f_in_12')
sm_f_in_13 = dae.add_z('sm_f_in_13')
sm_f_in_14 = dae.add_z('sm_f_in_14')
sm_f_in_15 = dae.add_z('sm_f_in_15')
sm_f_add_0 = dae.add_z('sm_f_add_0')
sm_f_add_2 = dae.add_z('sm_f_add_2')
sm_f_add_9 = dae.add_z('sm_f_add_9')
sm_f_add_14 = dae.add_z('sm_f_add_14')
sm_f_ex_0 = dae.add_z('sm_f_ex_0')
sm_f_ex_1 = dae.add_z('sm_f_ex_1')
sm_f_ex_2 = dae.add_z('sm_f_ex_2')
sm_f_ex_3 = dae.add_z('sm_f_ex_3')
sm_f_ex_4 = dae.add_z('sm_f_ex_4')
sm_f_ex_5 = dae.add_z('sm_f_ex_5')
sm_f_ex_6 = dae.add_z('sm_f_ex_6')
sm_f_ex_7 = dae.add_z('sm_f_ex_7')
sm_f_ex_8 = dae.add_z('sm_f_ex_8')
sm_f_ex_9 = dae.add_z('sm_f_ex_9')
sm_f_ex_10 = dae.add_z('sm_f_ex_10')
sm_f_ex_11 = dae.add_z('sm_f_ex_11')
sm_f_ex_12 = dae.add_z('sm_f_ex_12')
sm_f_ex_13 = dae.add_z('sm_f_ex_13')
sm_f_ex_14 = dae.add_z('sm_f_ex_14')
sm_f_ex_15 = dae.add_z('sm_f_ex_15')
sm_f_gs_1 = dae.add_z('sm_f_gs_1')
sm_f_gs_2 = dae.add_z('sm_f_gs_2')
sm_f_mm_0 = dae.add_z('sm_f_mm_0')
sm_f_mm_1 = dae.add_z('sm_f_mm_1')
sm_f_mm_2 = dae.add_z('sm_f_mm_2')
sm_f_mm_3 = dae.add_z('sm_f_mm_3')
sm_f_mm_4 = dae.add_z('sm_f_mm_4')
sm_f_mm_5 = dae.add_z('sm_f_mm_5')
sm_ntx = dae.add_z('sm_ntx')
sm_yx_0 = dae.add_z('sm_yx_0')
sm_yx_1 = dae.add_z('sm_yx_1')
sm_yx_2 = dae.add_z('sm_yx_2')
sm_yx_3 = dae.add_z('sm_yx_3')
sm_yx_4 = dae.add_z('sm_yx_4')
sm_yx_5 = dae.add_z('sm_yx_5')
sm_yx_6 = dae.add_z('sm_yx_6')
sm_yx_7 = dae.add_z('sm_yx_7')
sm_yx_8 = dae.add_z('sm_yx_8')
sm_yx_9 = dae.add_z('sm_yx_9')
sm_yx_10 = dae.add_z('sm_yx_10')
sm_yx_11 = dae.add_z('sm_yx_11')
sm_yx_12 = dae.add_z('sm_yx_12')
sm_y_prime_0 = dae.add_z('sm_y_prime_0')
sm_y_prime_1 = dae.add_z('sm_y_prime_1')
sm_y_prime_2 = dae.add_z('sm_y_prime_2')
sm_y_prime_3 = dae.add_z('sm_y_prime_3')
sm_y_prime_4 = dae.add_z('sm_y_prime_4')
sm_y_prime_5 = dae.add_z('sm_y_prime_5')
sm_y_prime_6 = dae.add_z('sm_y_prime_6')
sm_y_prime_7 = dae.add_z('sm_y_prime_7')
sm_y_prime_8 = dae.add_z('sm_y_prime_8')
sm_y_prime_9 = dae.add_z('sm_y_prime_9')
sm_y_prime_10 = dae.add_z('sm_y_prime_10')
sm_y_prime_11 = dae.add_z('sm_y_prime_11')
sm_y_prime_12 = dae.add_z('sm_y_prime_12')
sm_y_prime_13 = dae.add_z('sm_y_prime_13')
sm_y_prime_14 = dae.add_z('sm_y_prime_14')
sm_y_prime_15 = dae.add_z('sm_y_prime_15')
sm_m_0 = dae.add_z('sm_m_0')
sm_m_1 = dae.add_z('sm_m_1')
sm_m_2 = dae.add_z('sm_m_2')
sm_m_3 = dae.add_z('sm_m_3')
sm_m_4 = dae.add_z('sm_m_4')
sm_m_5 = dae.add_z('sm_m_5')
sm_m_6 = dae.add_z('sm_m_6')
sm_m_7 = dae.add_z('sm_m_7')
sm_m_8 = dae.add_z('sm_m_8')
sm_m_9 = dae.add_z('sm_m_9')
sm_m_10 = dae.add_z('sm_m_10')
sm_m_11 = dae.add_z('sm_m_11')
sm_m_12 = dae.add_z('sm_m_12')
sm_m_13 = dae.add_z('sm_m_13')
sm_m_14 = dae.add_z('sm_m_14')
sm_m_15 = dae.add_z('sm_m_15')
sm_ym_0 = dae.add_z('sm_ym_0')
sm_ym_1 = dae.add_z('sm_ym_1')
sm_ym_2 = dae.add_z('sm_ym_2')
sm_ym_3 = dae.add_z('sm_ym_3')
sm_ym_4 = dae.add_z('sm_ym_4')
sm_ym_5 = dae.add_z('sm_ym_5')
sm_ym_6 = dae.add_z('sm_ym_6')
sm_ym_7 = dae.add_z('sm_ym_7')
sm_ym_8 = dae.add_z('sm_ym_8')
sm_ym_9 = dae.add_z('sm_ym_9')
sm_ym_10 = dae.add_z('sm_ym_10')
sm_ym_11 = dae.add_z('sm_ym_11')
sm_ym_12 = dae.add_z('sm_ym_12')
sm_ym_13 = dae.add_z('sm_ym_13')
sm_ym_14 = dae.add_z('sm_ym_14')
sm_ym_15 = dae.add_z('sm_ym_15')
sm_m_cao_in = dae.add_z('sm_m_cao_in')
sm_m_dol_in = dae.add_z('sm_m_dol_in')
sm_f_in_o2_nm3_lance = dae.add_z('sm_f_in_o2_nm3_lance')
sm_f_in_kgm_c_tot = dae.add_z('sm_f_in_kgm_c_tot')
sm_o2_switch1 = dae.add_z('sm_o2_switch1')
sm_o2_switch2 = dae.add_z('sm_o2_switch2')
sm_o2_switch3 = dae.add_z('sm_o2_switch3')
sm_jetslag1 = dae.add_z('sm_jetslag1')
sm_jetslag2 = dae.add_z('sm_jetslag2')
sm_jetslag3 = dae.add_z('sm_jetslag3')
sm_f_in_o2_nm3 = dae.add_z('sm_f_in_o2_nm3')
sm_f_c_charge = dae.add_z('sm_f_c_charge')
sm_cpdt_0 = dae.add_z('sm_cpdt_0')
sm_cpdt_1 = dae.add_z('sm_cpdt_1')
sm_cpdt_2 = dae.add_z('sm_cpdt_2')
sm_cpdt_3 = dae.add_z('sm_cpdt_3')
sm_cpdt_4 = dae.add_z('sm_cpdt_4')
sm_cpdt_5 = dae.add_z('sm_cpdt_5')
sm_cpdt_6 = dae.add_z('sm_cpdt_6')
sm_cpdt_7 = dae.add_z('sm_cpdt_7')
sm_cpdt_8 = dae.add_z('sm_cpdt_8')
sm_cpdt_9 = dae.add_z('sm_cpdt_9')
sm_cpdt_10 = dae.add_z('sm_cpdt_10')
sm_cpdt_11 = dae.add_z('sm_cpdt_11')
sm_cpdt_12 = dae.add_z('sm_cpdt_12')
sm_cpdt_13 = dae.add_z('sm_cpdt_13')
sm_cpdt_14 = dae.add_z('sm_cpdt_14')
sm_cpdt_15 = dae.add_z('sm_cpdt_15')
sm_cptdt_0 = dae.add_z('sm_cptdt_0')
sm_cptdt_1 = dae.add_z('sm_cptdt_1')
sm_cptdt_2 = dae.add_z('sm_cptdt_2')
sm_cptdt_3 = dae.add_z('sm_cptdt_3')
sm_cptdt_4 = dae.add_z('sm_cptdt_4')
sm_cptdt_5 = dae.add_z('sm_cptdt_5')
sm_cptdt_6 = dae.add_z('sm_cptdt_6')
sm_cptdt_7 = dae.add_z('sm_cptdt_7')
sm_cptdt_8 = dae.add_z('sm_cptdt_8')
sm_cptdt_9 = dae.add_z('sm_cptdt_9')
sm_cptdt_10 = dae.add_z('sm_cptdt_10')
sm_cptdt_11 = dae.add_z('sm_cptdt_11')
sm_cptdt_12 = dae.add_z('sm_cptdt_12')
sm_cptdt_13 = dae.add_z('sm_cptdt_13')
sm_t_mm = dae.add_z('sm_t_mm')
sm_v_gs = dae.add_z('sm_v_gs')
sm_t_gs = dae.add_z('sm_t_gs')
sm_hf = dae.add_z('sm_hf')
sm_sigma_foam = dae.add_z('sm_sigma_foam')
sm_hs = dae.add_z('sm_hs')
sm_dens_slag = dae.add_z('sm_dens_slag')
sm_fi_foam = dae.add_z('sm_fi_foam')
sm_surf_tension = dae.add_z('sm_surf_tension')
sm_mu = dae.add_z('sm_mu')
sm_e1 = dae.add_z('sm_e1')
sm_e2 = dae.add_z('sm_e2')
sm_ef = dae.add_z('sm_ef')
sm_m_solid = dae.add_z('sm_m_solid')
sm_b_mu = dae.add_z('sm_b_mu')
sm_alpha_mu = dae.add_z('sm_alpha_mu')
sm_xg_star = dae.add_z('sm_xg_star')
sm_xm_star = dae.add_z('sm_xm_star')
sm_xa_star = dae.add_z('sm_xa_star')
sm_ym_prime_0 = dae.add_z('sm_ym_prime_0')
sm_ym_prime_1 = dae.add_z('sm_ym_prime_1')
sm_ym_prime_2 = dae.add_z('sm_ym_prime_2')
sm_ym_prime_3 = dae.add_z('sm_ym_prime_3')
sm_ym_prime_4 = dae.add_z('sm_ym_prime_4')
sm_ym_prime_5 = dae.add_z('sm_ym_prime_5')
sm_ym_prime_6 = dae.add_z('sm_ym_prime_6')
sm_ym_prime_7 = dae.add_z('sm_ym_prime_7')
sm_ym_prime_8 = dae.add_z('sm_ym_prime_8')
sm_ym_prime_9 = dae.add_z('sm_ym_prime_9')
sm_ym_prime_10 = dae.add_z('sm_ym_prime_10')
sm_ym_prime_11 = dae.add_z('sm_ym_prime_11')
sm_ym_prime_12 = dae.add_z('sm_ym_prime_12')
sm_ym_prime_13 = dae.add_z('sm_ym_prime_13')
sm_ym_prime_14 = dae.add_z('sm_ym_prime_14')
sm_ym_prime_15 = dae.add_z('sm_ym_prime_15')
sm_q = dae.add_z('sm_q')
sm_q_mm_sm = dae.add_z('sm_q_mm_sm')
sm_h_0 = dae.add_z('sm_h_0')
sm_h_1 = dae.add_z('sm_h_1')
sm_h_2 = dae.add_z('sm_h_2')
sm_h_3 = dae.add_z('sm_h_3')
sm_h_4 = dae.add_z('sm_h_4')
sm_h_5 = dae.add_z('sm_h_5')
sm_h_6 = dae.add_z('sm_h_6')
sm_h_7 = dae.add_z('sm_h_7')
sm_h_8 = dae.add_z('sm_h_8')
sm_h_9 = dae.add_z('sm_h_9')
sm_h_10 = dae.add_z('sm_h_10')
sm_h_11 = dae.add_z('sm_h_11')
sm_h_12 = dae.add_z('sm_h_12')
sm_h_13 = dae.add_z('sm_h_13')
sm_h_14 = dae.add_z('sm_h_14')
sm_h_15 = dae.add_z('sm_h_15')
sm_h_mm_0 = dae.add_z('sm_h_mm_0')
sm_h_mm_1 = dae.add_z('sm_h_mm_1')
sm_h_mm_2 = dae.add_z('sm_h_mm_2')
sm_h_mm_3 = dae.add_z('sm_h_mm_3')
sm_h_mm_4 = dae.add_z('sm_h_mm_4')
sm_h_mm_5 = dae.add_z('sm_h_mm_5')
sm_hflow_in = dae.add_z('sm_hflow_in')
sm_hflow_out = dae.add_z('sm_hflow_out')
sm_h_flow_in_0 = dae.add_z('sm_h_flow_in_0')
sm_h_flow_in_1 = dae.add_z('sm_h_flow_in_1')
sm_h_flow_in_2 = dae.add_z('sm_h_flow_in_2')
sm_h_flow_in_3 = dae.add_z('sm_h_flow_in_3')
sm_h_flow_in_4 = dae.add_z('sm_h_flow_in_4')
sm_h_flow_in_5 = dae.add_z('sm_h_flow_in_5')
sm_h_flow_in_6 = dae.add_z('sm_h_flow_in_6')
sm_h_flow_in_7 = dae.add_z('sm_h_flow_in_7')
sm_h_flow_in_8 = dae.add_z('sm_h_flow_in_8')
sm_h_flow_in_9 = dae.add_z('sm_h_flow_in_9')
sm_h_flow_in_10 = dae.add_z('sm_h_flow_in_10')
sm_h_flow_in_11 = dae.add_z('sm_h_flow_in_11')
sm_h_flow_in_12 = dae.add_z('sm_h_flow_in_12')
sm_h_flow_in_13 = dae.add_z('sm_h_flow_in_13')
sm_h_flow_in_14 = dae.add_z('sm_h_flow_in_14')
sm_h_flow_in_15 = dae.add_z('sm_h_flow_in_15')
sm_cpdt_mm_0 = dae.add_z('sm_cpdt_mm_0')
sm_cpdt_mm_1 = dae.add_z('sm_cpdt_mm_1')
sm_cpdt_mm_2 = dae.add_z('sm_cpdt_mm_2')
sm_cpdt_mm_3 = dae.add_z('sm_cpdt_mm_3')
sm_cpdt_mm_4 = dae.add_z('sm_cpdt_mm_4')
sm_cpdt_mm_5 = dae.add_z('sm_cpdt_mm_5')
ss_cpdt_1 = dae.add_z('ss_cpdt_1')
ss_n_ex = dae.add_z('ss_n_ex')
ss_m_liquid = dae.add_z('ss_m_liquid')
ss_n = dae.add_z('ss_n')
ss_n_in = dae.add_z('ss_n_in')
ss_q_chg = dae.add_z('ss_q_chg')
ss_q_gas = dae.add_z('ss_q_gas')
ss_q_mm_ss = dae.add_z('ss_q_mm_ss')
ss_q_power_ss = dae.add_z('ss_q_power_ss')
ss_q_vap = dae.add_z('ss_q_vap')
ss_t_mm = dae.add_z('ss_t_mm')
ss_soliddum5 = dae.add_z('ss_soliddum5')
ss_t_gs = dae.add_z('ss_t_gs')
ss_e_f = dae.add_z('ss_e_f')
ss_q_rad_ss = dae.add_z('ss_q_rad_ss')
ss_q_wall = dae.add_z('ss_q_wall')
mm_y_c = dae.add_z('mm_y_c')

# Parameter declarations
sm_e2_beta = dae.add_p('sm_e2_beta')
sm_teta1 = dae.add_p('sm_teta1')
rd_a_roof = dae.add_p('rd_a_roof')
rd_rho_bulk = dae.add_p('rd_rho_bulk')
sm_bias_o2_sm_star = dae.add_p('sm_bias_o2_sm_star')
gs_pi = dae.add_p('gs_pi')
gs_t_melt = dae.add_p('gs_t_melt')
gs_t_add = dae.add_p('gs_t_add')
sm_r = dae.add_p('sm_r')
rd_scrapcharge1st = dae.add_p('rd_scrapcharge1st')
gs_n_og = dae.add_p('gs_n_og')
sm_alpha1 = dae.add_p('sm_alpha1')
sm_cs_a = dae.add_p('sm_cs_a')
sm_k_c = dae.add_p('sm_k_c')
gs_xoil = dae.add_p('gs_xoil')
sm_e1_alpha = dae.add_p('sm_e1_alpha')
ss_k_dt = dae.add_p('ss_k_dt')
mm_alpha = dae.add_p('mm_alpha')
ss_k_dm = dae.add_p('ss_k_dm')
rd_tube = dae.add_p('rd_tube')
ss_k_cool = dae.add_p('ss_k_cool')
rd_r_roof = dae.add_p('rd_r_roof')
sm_beta1 = dae.add_p('sm_beta1')
sm_pi = dae.add_p('sm_pi')
rd_vtot = dae.add_p('rd_vtot')
gs_dhvap = dae.add_p('gs_dhvap')
ss_t_melt = dae.add_p('ss_t_melt')
rd_a_bath = dae.add_p('rd_a_bath')
ss_rho_bulk = dae.add_p('ss_rho_bulk')
sm_fi_beta = dae.add_p('sm_fi_beta')
sm_xdol = dae.add_p('sm_xdol')
rd_h_furnace_t = dae.add_p('rd_h_furnace_t')
rd_rho_tube = dae.add_p('rd_rho_tube')
ss_abg = dae.add_p('ss_abg')
gs_r = dae.add_p('gs_r')
ss_mw = dae.add_p('ss_mw')
sm_e1_beta = dae.add_p('sm_e1_beta')
gs_ttref = dae.add_p('gs_ttref')
mm_sub = dae.add_p('mm_sub')
ss_alpha = dae.add_p('ss_alpha')
rd_v_furnace = dae.add_p('rd_v_furnace')
sm_ym_p2o5 = dae.add_p('sm_ym_p2o5')
sm_scrapcharge1st = dae.add_p('sm_scrapcharge1st')
mm_y_star_c = dae.add_p('mm_y_star_c')
gs_bias_o2_gs_star = dae.add_p('gs_bias_o2_gs_star')
sm_tref = dae.add_p('sm_tref')
gs_f_h2o_star = dae.add_p('gs_f_h2o_star')
rd_h_gs = dae.add_p('rd_h_gs')
gs_dumgas2 = dae.add_p('gs_dumgas2')
gs_tref = dae.add_p('gs_tref')
sm_e2_alpha = dae.add_p('sm_e2_alpha')
rd_cp_wall = dae.add_p('rd_cp_wall')
rd_rho_heel = dae.add_p('rd_rho_heel')
gs_ttadd = dae.add_p('gs_ttadd')
sm_t_add = dae.add_p('sm_t_add')
mm_k_p = dae.add_p('mm_k_p')
mm_k_mcool = dae.add_p('mm_k_mcool')
mm_teta_l = dae.add_p('mm_teta_l')
ss_cp = dae.add_p('ss_cp')
rd_cp_roof = dae.add_p('rd_cp_roof')
mm_k_m = dae.add_p('mm_k_m')
gs_ea3 = dae.add_p('gs_ea3')
gs_ea2 = dae.add_p('gs_ea2')
gs_ea1 = dae.add_p('gs_ea1')
sm_k_dc = dae.add_p('sm_k_dc')
rd_a_solid = dae.add_p('rd_a_solid')
sm_cao_fr_dol = dae.add_p('sm_cao_fr_dol')
rd_heel = dae.add_p('rd_heel')
gs_p_std = dae.add_p('gs_p_std')
mm_gamma_d = dae.add_p('mm_gamma_d')
mm_k_t2 = dae.add_p('mm_k_t2')
rd_ua_t1 = dae.add_p('rd_ua_t1')
ss_dh_fus = dae.add_p('ss_dh_fus')
gs_k_oil = dae.add_p('gs_k_oil')
gs_k_po2 = dae.add_p('gs_k_po2')
sm_p_std = dae.add_p('sm_p_std')
sm_k_cao = dae.add_p('sm_k_cao')
sm_furnace_dia = dae.add_p('sm_furnace_dia')
ss_k_t1 = dae.add_p('ss_k_t1')
ss_k_t3 = dae.add_p('ss_k_t3')
sm_xc_imp = dae.add_p('sm_xc_imp')
sm_teta_l = dae.add_p('sm_teta_l')
rd_k_steel = dae.add_p('rd_k_steel')
ss_k_p = dae.add_p('ss_k_p')
rd_pi = dae.add_p('rd_pi')
sm_xcao = dae.add_p('sm_xcao')
sm_fi_alpha = dae.add_p('sm_fi_alpha')
sm_cp_0 = dae.add_p('sm_cp_0')
sm_cp_1 = dae.add_p('sm_cp_1')
sm_cp_2 = dae.add_p('sm_cp_2')
sm_cp_3 = dae.add_p('sm_cp_3')
sm_cp_4 = dae.add_p('sm_cp_4')
sm_cp_5 = dae.add_p('sm_cp_5')
sm_cp_6 = dae.add_p('sm_cp_6')
sm_cp_7 = dae.add_p('sm_cp_7')
sm_cp_8 = dae.add_p('sm_cp_8')
sm_cp_9 = dae.add_p('sm_cp_9')
sm_cp_10 = dae.add_p('sm_cp_10')
sm_cp_11 = dae.add_p('sm_cp_11')
sm_cp_12 = dae.add_p('sm_cp_12')
sm_cp_13 = dae.add_p('sm_cp_13')
sm_cp_14 = dae.add_p('sm_cp_14')
sm_cp_15 = dae.add_p('sm_cp_15')
rd_t_ec_0 = dae.add_p('rd_t_ec_0')
rd_t_ec_1 = dae.add_p('rd_t_ec_1')
sm_mw_0 = dae.add_p('sm_mw_0')
sm_mw_1 = dae.add_p('sm_mw_1')
sm_mw_2 = dae.add_p('sm_mw_2')
sm_mw_3 = dae.add_p('sm_mw_3')
sm_mw_4 = dae.add_p('sm_mw_4')
sm_mw_5 = dae.add_p('sm_mw_5')
sm_mw_6 = dae.add_p('sm_mw_6')
sm_mw_7 = dae.add_p('sm_mw_7')
sm_mw_8 = dae.add_p('sm_mw_8')
sm_mw_9 = dae.add_p('sm_mw_9')
sm_mw_10 = dae.add_p('sm_mw_10')
sm_mw_11 = dae.add_p('sm_mw_11')
sm_mw_12 = dae.add_p('sm_mw_12')
sm_mw_13 = dae.add_p('sm_mw_13')
sm_mw_14 = dae.add_p('sm_mw_14')
sm_mw_15 = dae.add_p('sm_mw_15')
mm_x_0 = dae.add_p('mm_x_0')
mm_x_1 = dae.add_p('mm_x_1')
mm_x_2 = dae.add_p('mm_x_2')
mm_x_3 = dae.add_p('mm_x_3')
mm_x_4 = dae.add_p('mm_x_4')
mm_x_5 = dae.add_p('mm_x_5')
sm_v_molar_0 = dae.add_p('sm_v_molar_0')
sm_v_molar_1 = dae.add_p('sm_v_molar_1')
sm_v_molar_2 = dae.add_p('sm_v_molar_2')
sm_v_molar_3 = dae.add_p('sm_v_molar_3')
sm_v_molar_4 = dae.add_p('sm_v_molar_4')
sm_v_molar_5 = dae.add_p('sm_v_molar_5')
sm_v_molar_6 = dae.add_p('sm_v_molar_6')
sm_v_molar_7 = dae.add_p('sm_v_molar_7')
sm_v_molar_8 = dae.add_p('sm_v_molar_8')
sm_v_molar_9 = dae.add_p('sm_v_molar_9')
sm_v_molar_10 = dae.add_p('sm_v_molar_10')
sm_v_molar_11 = dae.add_p('sm_v_molar_11')
sm_v_molar_12 = dae.add_p('sm_v_molar_12')
sm_v_molar_13 = dae.add_p('sm_v_molar_13')
sm_v_molar_14 = dae.add_p('sm_v_molar_14')
sm_v_molar_15 = dae.add_p('sm_v_molar_15')
sm_cpdt_add_0 = dae.add_p('sm_cpdt_add_0')
sm_cpdt_add_1 = dae.add_p('sm_cpdt_add_1')
sm_cpdt_add_2 = dae.add_p('sm_cpdt_add_2')
sm_cpdt_add_3 = dae.add_p('sm_cpdt_add_3')
sm_cpdt_add_4 = dae.add_p('sm_cpdt_add_4')
sm_cpdt_add_5 = dae.add_p('sm_cpdt_add_5')
sm_cpdt_add_6 = dae.add_p('sm_cpdt_add_6')
sm_cpdt_add_7 = dae.add_p('sm_cpdt_add_7')
sm_cpdt_add_8 = dae.add_p('sm_cpdt_add_8')
sm_cpdt_add_9 = dae.add_p('sm_cpdt_add_9')
sm_cpdt_add_10 = dae.add_p('sm_cpdt_add_10')
sm_cpdt_add_11 = dae.add_p('sm_cpdt_add_11')
sm_cpdt_add_12 = dae.add_p('sm_cpdt_add_12')
sm_cpdt_add_13 = dae.add_p('sm_cpdt_add_13')
sm_cpdt_add_14 = dae.add_p('sm_cpdt_add_14')
sm_cpdt_add_15 = dae.add_p('sm_cpdt_add_15')
gs_mw_0 = dae.add_p('gs_mw_0')
gs_mw_1 = dae.add_p('gs_mw_1')
gs_mw_2 = dae.add_p('gs_mw_2')
gs_mw_3 = dae.add_p('gs_mw_3')
gs_mw_4 = dae.add_p('gs_mw_4')
gs_mw_5 = dae.add_p('gs_mw_5')
gs_mw_6 = dae.add_p('gs_mw_6')
gs_mw_7 = dae.add_p('gs_mw_7')
gs_x_air_0 = dae.add_p('gs_x_air_0')
gs_x_air_1 = dae.add_p('gs_x_air_1')
gs_x_air_2 = dae.add_p('gs_x_air_2')
gs_x_air_3 = dae.add_p('gs_x_air_3')
gs_x_air_4 = dae.add_p('gs_x_air_4')
gs_x_air_5 = dae.add_p('gs_x_air_5')
gs_x_air_6 = dae.add_p('gs_x_air_6')
gs_x_air_7 = dae.add_p('gs_x_air_7')
gs_h_add_0 = dae.add_p('gs_h_add_0')
gs_h_add_1 = dae.add_p('gs_h_add_1')
gs_h_add_2 = dae.add_p('gs_h_add_2')
gs_h_add_3 = dae.add_p('gs_h_add_3')
gs_h_add_4 = dae.add_p('gs_h_add_4')
gs_h_add_5 = dae.add_p('gs_h_add_5')
gs_h_add_6 = dae.add_p('gs_h_add_6')
gs_h_add_7 = dae.add_p('gs_h_add_7')
mm_beta_i_0 = dae.add_p('mm_beta_i_0')
mm_beta_i_1 = dae.add_p('mm_beta_i_1')
mm_beta_i_2 = dae.add_p('mm_beta_i_2')
mm_beta_i_3 = dae.add_p('mm_beta_i_3')
mm_beta_i_4 = dae.add_p('mm_beta_i_4')
mm_beta_i_5 = dae.add_p('mm_beta_i_5')
gs_cp_0 = dae.add_p('gs_cp_0')
gs_cp_1 = dae.add_p('gs_cp_1')
gs_cp_2 = dae.add_p('gs_cp_2')
gs_cp_3 = dae.add_p('gs_cp_3')
gs_cp_4 = dae.add_p('gs_cp_4')
gs_cp_5 = dae.add_p('gs_cp_5')
gs_cp_6 = dae.add_p('gs_cp_6')
gs_cp_7 = dae.add_p('gs_cp_7')
mm_abc_0 = dae.add_p('mm_abc_0')
mm_abc_1 = dae.add_p('mm_abc_1')
mm_abc_2 = dae.add_p('mm_abc_2')
mm_abc_3 = dae.add_p('mm_abc_3')
mm_abc_4 = dae.add_p('mm_abc_4')
mm_abc_5 = dae.add_p('mm_abc_5')
mm_abg_0 = dae.add_p('mm_abg_0')
mm_abg_1 = dae.add_p('mm_abg_1')
mm_abg_2 = dae.add_p('mm_abg_2')
mm_abg_3 = dae.add_p('mm_abg_3')
mm_abg_4 = dae.add_p('mm_abg_4')
mm_abg_5 = dae.add_p('mm_abg_5')
rd_t_xc_0 = dae.add_p('rd_t_xc_0')
rd_t_xc_1 = dae.add_p('rd_t_xc_1')
sm_dg0_0 = dae.add_p('sm_dg0_0')
sm_dg0_1 = dae.add_p('sm_dg0_1')
sm_dg0_2 = dae.add_p('sm_dg0_2')
sm_dg0_3 = dae.add_p('sm_dg0_3')
sm_dg0_4 = dae.add_p('sm_dg0_4')
sm_dg0_5 = dae.add_p('sm_dg0_5')
sm_dg0_6 = dae.add_p('sm_dg0_6')
sm_dg0_7 = dae.add_p('sm_dg0_7')
sm_dg0_8 = dae.add_p('sm_dg0_8')
sm_dg0_9 = dae.add_p('sm_dg0_9')
sm_dg0_10 = dae.add_p('sm_dg0_10')
sm_dg0_11 = dae.add_p('sm_dg0_11')
sm_dg0_12 = dae.add_p('sm_dg0_12')
sm_dg0_13 = dae.add_p('sm_dg0_13')
sm_dg0_14 = dae.add_p('sm_dg0_14')
sm_dg0_15 = dae.add_p('sm_dg0_15')
mm_mw_0 = dae.add_p('mm_mw_0')
mm_mw_1 = dae.add_p('mm_mw_1')
mm_mw_2 = dae.add_p('mm_mw_2')
mm_mw_3 = dae.add_p('mm_mw_3')
mm_mw_4 = dae.add_p('mm_mw_4')
mm_mw_5 = dae.add_p('mm_mw_5')
mm_cp_0 = dae.add_p('mm_cp_0')
mm_cp_1 = dae.add_p('mm_cp_1')
mm_cp_2 = dae.add_p('mm_cp_2')
mm_cp_3 = dae.add_p('mm_cp_3')
mm_cp_4 = dae.add_p('mm_cp_4')
mm_cp_5 = dae.add_p('mm_cp_5')
sm_dhf0_0 = dae.add_p('sm_dhf0_0')
sm_dhf0_1 = dae.add_p('sm_dhf0_1')
sm_dhf0_2 = dae.add_p('sm_dhf0_2')
sm_dhf0_3 = dae.add_p('sm_dhf0_3')
sm_dhf0_4 = dae.add_p('sm_dhf0_4')
sm_dhf0_5 = dae.add_p('sm_dhf0_5')
sm_dhf0_6 = dae.add_p('sm_dhf0_6')
sm_dhf0_7 = dae.add_p('sm_dhf0_7')
sm_dhf0_8 = dae.add_p('sm_dhf0_8')
sm_dhf0_9 = dae.add_p('sm_dhf0_9')
sm_dhf0_10 = dae.add_p('sm_dhf0_10')
sm_dhf0_11 = dae.add_p('sm_dhf0_11')
sm_dhf0_12 = dae.add_p('sm_dhf0_12')
sm_dhf0_13 = dae.add_p('sm_dhf0_13')
sm_dhf0_14 = dae.add_p('sm_dhf0_14')
sm_dhf0_15 = dae.add_p('sm_dhf0_15')
gs_dhf0_0 = dae.add_p('gs_dhf0_0')
gs_dhf0_1 = dae.add_p('gs_dhf0_1')
gs_dhf0_2 = dae.add_p('gs_dhf0_2')
gs_dhf0_3 = dae.add_p('gs_dhf0_3')
gs_dhf0_4 = dae.add_p('gs_dhf0_4')
gs_dhf0_5 = dae.add_p('gs_dhf0_5')
gs_dhf0_6 = dae.add_p('gs_dhf0_6')
gs_dhf0_7 = dae.add_p('gs_dhf0_7')
sm_dens_i_0 = dae.add_p('sm_dens_i_0')
sm_dens_i_1 = dae.add_p('sm_dens_i_1')
sm_dens_i_2 = dae.add_p('sm_dens_i_2')
sm_dens_i_3 = dae.add_p('sm_dens_i_3')
sm_dens_i_4 = dae.add_p('sm_dens_i_4')
sm_dens_i_5 = dae.add_p('sm_dens_i_5')
sm_dens_i_6 = dae.add_p('sm_dens_i_6')
sm_dens_i_7 = dae.add_p('sm_dens_i_7')
sm_dens_i_8 = dae.add_p('sm_dens_i_8')
sm_dens_i_9 = dae.add_p('sm_dens_i_9')
sm_dens_i_10 = dae.add_p('sm_dens_i_10')
sm_dens_i_11 = dae.add_p('sm_dens_i_11')
sm_dens_i_12 = dae.add_p('sm_dens_i_12')
sm_dens_i_13 = dae.add_p('sm_dens_i_13')
sm_dens_i_14 = dae.add_p('sm_dens_i_14')
sm_dens_i_15 = dae.add_p('sm_dens_i_15')
gs_dg0_0 = dae.add_p('gs_dg0_0')
gs_dg0_1 = dae.add_p('gs_dg0_1')
gs_dg0_2 = dae.add_p('gs_dg0_2')
gs_dg0_3 = dae.add_p('gs_dg0_3')
gs_dg0_4 = dae.add_p('gs_dg0_4')
gs_dg0_5 = dae.add_p('gs_dg0_5')
gs_dg0_6 = dae.add_p('gs_dg0_6')
gs_dg0_7 = dae.add_p('gs_dg0_7')
sm_h_add_0 = dae.add_p('sm_h_add_0')
sm_h_add_1 = dae.add_p('sm_h_add_1')
sm_h_add_2 = dae.add_p('sm_h_add_2')
sm_h_add_3 = dae.add_p('sm_h_add_3')
sm_h_add_4 = dae.add_p('sm_h_add_4')
sm_h_add_5 = dae.add_p('sm_h_add_5')
sm_h_add_6 = dae.add_p('sm_h_add_6')
sm_h_add_7 = dae.add_p('sm_h_add_7')
sm_h_add_8 = dae.add_p('sm_h_add_8')
sm_h_add_9 = dae.add_p('sm_h_add_9')
sm_h_add_10 = dae.add_p('sm_h_add_10')
sm_h_add_11 = dae.add_p('sm_h_add_11')
sm_h_add_12 = dae.add_p('sm_h_add_12')
sm_h_add_13 = dae.add_p('sm_h_add_13')
sm_h_add_14 = dae.add_p('sm_h_add_14')
sm_h_add_15 = dae.add_p('sm_h_add_15')
gs_cpdt_add_0 = dae.add_p('gs_cpdt_add_0')
gs_cpdt_add_1 = dae.add_p('gs_cpdt_add_1')
gs_cpdt_add_2 = dae.add_p('gs_cpdt_add_2')
gs_cpdt_add_3 = dae.add_p('gs_cpdt_add_3')
gs_cpdt_add_4 = dae.add_p('gs_cpdt_add_4')
gs_cpdt_add_5 = dae.add_p('gs_cpdt_add_5')
gs_cpdt_add_6 = dae.add_p('gs_cpdt_add_6')
gs_cpdt_add_7 = dae.add_p('gs_cpdt_add_7')
gs_n_sm_1 = dae.add_p('gs_n_sm_1')
gs_n_sm_3 = dae.add_p('gs_n_sm_3')
gs_n_sm_4 = dae.add_p('gs_n_sm_4')
gs_n_sm_5 = dae.add_p('gs_n_sm_5')
gs_n_sm_6 = dae.add_p('gs_n_sm_6')
gs_n_sm_7 = dae.add_p('gs_n_sm_7')
sm_f_add_1 = dae.add_p('sm_f_add_1')
sm_f_gs_0 = dae.add_p('sm_f_gs_0')
sm_f_gs_3 = dae.add_p('sm_f_gs_3')
sm_f_gs_4 = dae.add_p('sm_f_gs_4')
sm_f_gs_5 = dae.add_p('sm_f_gs_5')
sm_f_gs_6 = dae.add_p('sm_f_gs_6')
sm_f_gs_7 = dae.add_p('sm_f_gs_7')
sm_f_gs_8 = dae.add_p('sm_f_gs_8')
sm_f_gs_9 = dae.add_p('sm_f_gs_9')
sm_f_gs_10 = dae.add_p('sm_f_gs_10')
sm_f_gs_11 = dae.add_p('sm_f_gs_11')
sm_f_gs_12 = dae.add_p('sm_f_gs_12')
sm_f_gs_13 = dae.add_p('sm_f_gs_13')
sm_f_gs_14 = dae.add_p('sm_f_gs_14')
sm_f_gs_15 = dae.add_p('sm_f_gs_15')
sm_f_add_3 = dae.add_p('sm_f_add_3')
sm_f_add_4 = dae.add_p('sm_f_add_4')
sm_f_add_5 = dae.add_p('sm_f_add_5')
sm_f_add_6 = dae.add_p('sm_f_add_6')
sm_f_add_7 = dae.add_p('sm_f_add_7')
sm_f_add_8 = dae.add_p('sm_f_add_8')
sm_f_add_10 = dae.add_p('sm_f_add_10')
sm_f_add_11 = dae.add_p('sm_f_add_11')
sm_f_add_12 = dae.add_p('sm_f_add_12')
sm_f_add_13 = dae.add_p('sm_f_add_13')
sm_f_add_15 = dae.add_p('sm_f_add_15')
gs_f_in_nm3_0 = dae.add_p('gs_f_in_nm3_0')
gs_f_in_nm3_1 = dae.add_p('gs_f_in_nm3_1')
gs_f_in_nm3_4 = dae.add_p('gs_f_in_nm3_4')
gs_f_in_nm3_7 = dae.add_p('gs_f_in_nm3_7')
sm_b_8 = dae.add_p('sm_b_8')
mm_n_5 = dae.add_p('mm_n_5')
gs_n_ex_6 = dae.add_p('gs_n_ex_6')

# Input variable declarations (which are part of the parameters)
gs_f_in_kgm_oil         = dae.add_u('gs_f_in_kgm_oil')
sm_f_in_kgm_cao         = dae.add_u('sm_f_in_kgm_cao')
sm_f_in_kgm_c_charge    = dae.add_u('sm_f_in_kgm_c_charge')
sm_f_in_kgm_c_lance     = dae.add_u('sm_f_in_kgm_c_lance')
gs_f_h20_electrode      = dae.add_u('gs_f_h20_electrode')
ss_m_in_t               = dae.add_u('ss_m_in_t')
sm_f_in_kgm_dolo        = dae.add_u('sm_f_in_kgm_dolo')
gs_f_in_ch4             = dae.add_u('gs_f_in_ch4')
ss_p_arc                = dae.add_u('ss_p_arc')
sm_f_in_o2_jetbox2      = dae.add_u('sm_f_in_o2_jetbox2')
sm_f_in_o2_jetbox3      = dae.add_u('sm_f_in_o2_jetbox3')
sm_f_in_o2_jetbox1      = dae.add_u('sm_f_in_o2_jetbox1')

# Model noise declarations
w_0 = dae.add_u('w_0')
w_1 = dae.add_u('w_1')
w_2 = dae.add_u('w_2')
w_3 = dae.add_u('w_3')
w_4 = dae.add_u('w_4')
w_5 = dae.add_u('w_5')
w_6 = dae.add_u('w_6')
w_7 = dae.add_u('w_7')
w_8 = dae.add_u('w_8')
w_9 = dae.add_u('w_9')
w_10 = dae.add_u('w_10')
w_11 = dae.add_u('w_11')
w_12 = dae.add_u('w_12')
w_13 = dae.add_u('w_13')
w_14 = dae.add_u('w_14')
w_15 = dae.add_u('w_15')
w_16 = dae.add_u('w_16')
w_17 = dae.add_u('w_17')
w_18 = dae.add_u('w_18')
w_19 = dae.add_u('w_19')
w_20 = dae.add_u('w_20')
w_21 = dae.add_u('w_21')
w_22 = dae.add_u('w_22')
w_23 = dae.add_u('w_23')
w_24 = dae.add_u('w_24')
w_25 = dae.add_u('w_25')
w_26 = dae.add_u('w_26')
w_27 = dae.add_u('w_27')

# Add the model noise for the disturbance state
if disturb_add == 1:
    w_28 = dae.add_u('w_28')

# Scaling for sm_b_3
sm_b_3_scale = 100

# Scaling for sm_e
sm_e_scale = 0.001
               
# Declare differential equations in explicit form
if disturb_add ==1:
    ode = [(sm_f_in_0-sm_f_ex_1) + w_0,
            (((2*(sm_f_in_2-sm_f_ex_2))-sm_f_ex_1)+sm_f_in_9) + w_1,#
            sm_f_in_3 + w_2,
            #sm_f_in_6 + w_3,
            (sm_f_in_6)*sm_b_3_scale + w_3 + d_sm_b_3, # distubance state added
            (sm_f_in_8+sm_f_in_9) + w_4,#
            sm_f_in_10 + w_5,
            sm_f_in_12 + w_6,
            sm_f_in_14 + w_7,#
            ((sm_xcao*sm_f_in_kgm_cao)-(sm_k_cao*sm_m_cao)) + w_8,
            (((1-sm_xc_imp)*sm_f_in_kgm_c_charge)-(sm_k_dc*sm_m_c_float)) + w_9,
            ((sm_xdol*sm_f_in_kgm_dolo)-(sm_k_cao*sm_m_dol)) + w_10,
            ((sm_q+sm_hflow_in)-sm_hflow_out)*sm_e_scale + w_11,#
            (ss_m_in_t-((ss_n_ex*ss_mw)/1000)) + 0.00001*w_12,
            (((60*(((((ss_q_rad_ss+ss_q_power_ss)+ss_q_mm_ss)+ss_q_gas)-ss_q_chg)-ss_q_vap))*(1-(ss_t/ss_t_melt)))/((ss_cp*ss_n)*ss_k_dt+1e-7)) + 0.00001*w_13,#+0.0000001)) + w_13,
            ((mm_melt_0+mm_f_in)-mm_f_ex_0) + w_14,
            (mm_melt_1-mm_f_ex_1) + w_15,
            (mm_melt_2-mm_f_ex_2) + w_16,
            (mm_melt_3-mm_f_ex_3) + w_17,
            (mm_melt_4-mm_f_ex_4) + w_18,
            ((60*((((mm_q_rad+mm_p_power)-mm_q_mm_ss)-mm_q_cool)-mm_q_mm_sm))/(mm_sub*((((((mm_n_0*mm_cp_0)+(mm_n_1*mm_cp_1))+(mm_n_2*mm_cp_2))+(mm_n_3*mm_cp_3))+(mm_n_4*mm_cp_4))+(mm_n_5*mm_cp_5))+0.0000001)) + w_19,
            ((((gs_n_in_0-gs_n_ex_0)-gs_n_ex_1)+(gs_n_in_3-gs_n_ex_3))+(9*gs_n_in_6)) + w_20,#
            ((((gs_n_in_0-gs_n_ex_0)-(2*gs_n_ex_1))+(2*(gs_n_in_2-gs_n_ex_2)))+(gs_n_in_5-gs_n_ex_5)) + w_21,#
            ((((4*(gs_n_in_3-gs_n_ex_3))-(2*gs_n_ex_4))+(2*(gs_n_in_5-gs_n_ex_5)))+(20*gs_n_in_6)) + w_22,
            (2*(gs_n_in_7-gs_n_ex_7)) + w_23,
            ((gs_xoil*gs_n_oil_in)-gs_n_add_6) + w_24,
            (((gs_q+gs_fh_in)-gs_fh_out)/(0.03*gs_ntx+0.00000001)) + w_25,#
            (60*((((rd_rho_tube*rd_tube)*(rd_a_roof+0.01))*rd_cp_roof)*((rd_q_radroof+rd_q_gascool_0)-rd_q_cw_0))) + w_26,
            (60*((((rd_rho_tube*rd_tube)*(rd_a_wall+0.01))*rd_cp_wall)*((rd_q_gascool_1-rd_q_cw_1)+rd_q_radwall))) + w_27,
            (0 + w_28) # disturbace state equation
    ]
else:    
    ode = [(sm_f_in_0-sm_f_ex_1) + w_0,
            (((2*(sm_f_in_2-sm_f_ex_2))-sm_f_ex_1)+sm_f_in_9) + w_1,#
            sm_f_in_3 + w_2,
            #sm_f_in_6 + w_3,
            (sm_f_in_6)*sm_b_3_scale + w_3,# + d_sm_b_3, # distubance state added
            (sm_f_in_8+sm_f_in_9) + w_4,#
            sm_f_in_10 + w_5,
            sm_f_in_12 + w_6,
            sm_f_in_14 + w_7,#
            ((sm_xcao*sm_f_in_kgm_cao)-(sm_k_cao*sm_m_cao)) + w_8,
            (((1-sm_xc_imp)*sm_f_in_kgm_c_charge)-(sm_k_dc*sm_m_c_float)) + w_9,
            ((sm_xdol*sm_f_in_kgm_dolo)-(sm_k_cao*sm_m_dol)) + w_10,
            ((sm_q+sm_hflow_in)-sm_hflow_out)*sm_e_scale + w_11,#
            (ss_m_in_t-((ss_n_ex*ss_mw)/1000)) + w_12,
            (((60*(((((ss_q_rad_ss+ss_q_power_ss)+ss_q_mm_ss)+ss_q_gas)-ss_q_chg)-ss_q_vap))*(1-(ss_t/ss_t_melt)))/((ss_cp*ss_n)*ss_k_dt+0.0000001)) + w_13,
            ((mm_melt_0+mm_f_in)-mm_f_ex_0) + w_14,
            (mm_melt_1-mm_f_ex_1) + w_15,
            (mm_melt_2-mm_f_ex_2) + w_16,
            (mm_melt_3-mm_f_ex_3) + w_17,
            (mm_melt_4-mm_f_ex_4) + w_18,
            ((60*((((mm_q_rad+mm_p_power)-mm_q_mm_ss)-mm_q_cool)-mm_q_mm_sm))/(mm_sub*((((((mm_n_0*mm_cp_0)+(mm_n_1*mm_cp_1))+(mm_n_2*mm_cp_2))+(mm_n_3*mm_cp_3))+(mm_n_4*mm_cp_4))+(mm_n_5*mm_cp_5))+0.0000001)) + w_19,
            ((((gs_n_in_0-gs_n_ex_0)-gs_n_ex_1)+(gs_n_in_3-gs_n_ex_3))+(9*gs_n_in_6)) + w_20,#
            ((((gs_n_in_0-gs_n_ex_0)-(2*gs_n_ex_1))+(2*(gs_n_in_2-gs_n_ex_2)))+(gs_n_in_5-gs_n_ex_5)) + w_21,#
            ((((4*(gs_n_in_3-gs_n_ex_3))-(2*gs_n_ex_4))+(2*(gs_n_in_5-gs_n_ex_5)))+(20*gs_n_in_6)) + w_22,
            (2*(gs_n_in_7-gs_n_ex_7)) + w_23,
            ((gs_xoil*gs_n_oil_in)-gs_n_add_6) + w_24,
            (((gs_q+gs_fh_in)-gs_fh_out)/(0.03*gs_ntx+0.00000001)) + w_25,#
            (60*((((rd_rho_tube*rd_tube)*(rd_a_roof+0.01))*rd_cp_roof)*((rd_q_radroof+rd_q_gascool_0)-rd_q_cw_0))) + w_26,
            (60*((((rd_rho_tube*rd_tube)*(rd_a_wall+0.01))*rd_cp_wall)*((rd_q_gascool_1-rd_q_cw_1)+rd_q_radwall))) + w_27,
            #(0 + w_28) # disturbace state equation
    ]

# Declare algebraic equations in explicit form
# Had to divide the contents into 2 parts because python not taking more than
# 255 arguments for a function              
alg =[(mm_t-ss_t_mm),
          (mm_p_arc-ss_p_arc),
          (mm_q_mm_ss-ss_q_mm_ss),
          ((((((mm_m_mm_0+mm_m_mm_1)+mm_m_mm_2)+mm_m_mm_3)+mm_m_mm_4)+mm_m_mm_5)-ss_m_liquid),
          (sm_ef-ss_e_f),
          (((gs_dhvap*gs_n_add_6)/60)-ss_q_vap),
          ((((ss_m_ss*ss_k_t3)*(gs_n_add_2+gs_n_add_3))*(ss_t_gs-ss_t))-ss_q_gas),
          (gs_t-ss_t_gs),
          ((((((((((((((((sm_m_0+sm_m_1)+sm_m_2)+sm_m_3)+sm_m_4)+sm_m_5)+sm_m_6)+sm_m_7)+sm_m_8)+sm_m_9)+sm_m_10)+sm_m_11)+sm_m_12)+sm_m_13)+sm_m_14)+sm_m_15)-mm_m_slag),
          (sm_t-mm_t_sm),
          (sm_f_in_kgm_c_lance-mm_f_in_kgm_c_lance),
          (ss_n_ex-mm_f_melt_ss),
          (ss_m_ss-mm_m_solid),
          (sm_ef-mm_e_f),
          ((sm_f_in_o2_nm3/60)-mm_f_o2_lnc),
          (mm_y_slag_0-sm_yx_0),
          (mm_y_slag_1-sm_yx_3),
          (mm_y_slag_2-sm_yx_6),
          (mm_y_slag_3-sm_yx_10),
          (mm_y_slag_4-sm_yx_12),
          (mm_y_slag_5-sm_yx_8),
          (ss_q_wall-mm_q_wall),
          (mm_f_sm_0-sm_f_mm_0),
          (mm_f_sm_1-sm_f_mm_1),
          (mm_f_sm_2-sm_f_mm_2),
          (mm_f_sm_3-sm_f_mm_3),
          (mm_f_sm_4-sm_f_mm_4),
          (mm_f_sm_5-sm_f_mm_5),
          (gs_n_sm_2-sm_f_gs_2),
          (mm_t-sm_t_mm),
          (gs_t-sm_t_gs),
          (mm_m_solid-sm_m_solid),
          (mm_q_mm_sm-sm_q_mm_sm),
          (ss_t-gs_t_ss),
          (sm_f_in_o2_jetbox1-gs_f_in_o2_jetbox1),
          (sm_f_in_o2_jetbox2-gs_f_in_o2_jetbox2),
          (sm_f_in_o2_jetbox3-gs_f_in_o2_jetbox3),
          (sm_o2_switch1-gs_o2_switch1),
          (sm_o2_switch2-gs_o2_switch2),
          (sm_o2_switch3-gs_o2_switch3),
          (sm_f_gs_1-gs_n_sm_0),
          (((exp(sm_xl_2)-exp(gs_xl_2))*gs_k_po2)-(1000*gs_n_sm_2)),
          (sm_t-gs_t_sm),
          (ss_q_gas-gs_q_ss),
          ((rd_q_gascool_0+rd_q_gascool_1)-gs_q_cool),
          (gs_t-rd_t_gs),
          (mm_m_solid-rd_m_solid),
          (ss_m_liquid-rd_m_liquid),
          (ss_q_wall-rd_q_wall),
          (ss_e_f-rd_e_f),
          (sm_lambda_c-sm_sum_lambda_0),
          ((sm_lambda_c-sm_sum_lambda_1)+sm_lambda_o),
          ((2*sm_lambda_o)-sm_sum_lambda_2),
          (sm_lambda_fe-sm_sum_lambda_3),
          ((sm_lambda_fe-sm_sum_lambda_4)+sm_lambda_o),
          (((2*sm_lambda_fe)-sm_sum_lambda_5)+(3*sm_lambda_o)),
          (sm_lambda_mn-sm_sum_lambda_6),
          ((sm_lambda_mn-sm_sum_lambda_7)+sm_lambda_o),
          (sm_lambda_mg-sm_sum_lambda_8),
          ((sm_lambda_mg-sm_sum_lambda_9)+sm_lambda_o),
          (sm_lambda_si-sm_sum_lambda_10),
          ((sm_lambda_si-sm_sum_lambda_11)+(2*sm_lambda_o)),
          (sm_lambda_al-sm_sum_lambda_12),
          (((2*sm_lambda_al)+(3*sm_lambda_o))-sm_sum_lambda_13),
          ((sm_dg_0+((sm_r*sm_t)*sm_xl_0))+sm_sum_lambda_0),
          (((sm_dg_1+((sm_r*sm_t)*sm_xl_1))+sm_sum_lambda_1)),
          ((sm_dg_2+((sm_r*sm_t)*sm_xl_2))+sm_sum_lambda_2),
          ((sm_dg_3+((sm_r*sm_t)*sm_xl_3))+sm_sum_lambda_3),
          ((sm_dg_4+((sm_r*sm_t)*sm_xl_4))+sm_sum_lambda_4),
          ((sm_dg_5+((sm_r*sm_t)*sm_xl_5))+sm_sum_lambda_5),
          ((sm_dg_6+((sm_r*sm_t)*sm_xl_6))+sm_sum_lambda_6),
          ((sm_dg_7+((sm_r*sm_t)*sm_xl_7))+sm_sum_lambda_7),
          ((sm_dg_8+((sm_r*sm_t)*sm_xl_8))+sm_sum_lambda_8),
          ((sm_dg_9+((sm_r*sm_t)*sm_xl_9))+sm_sum_lambda_9),
          ((sm_dg_10+((sm_r*sm_t)*sm_xl_10))+sm_sum_lambda_10),
          ((sm_dg_11+((sm_r*sm_t)*sm_xl_11))+sm_sum_lambda_11),
          ((sm_dg_12+((sm_r*sm_t)*sm_xl_12))+sm_sum_lambda_12),
          ((sm_dg_13+((sm_r*sm_t)*sm_xl_13))+sm_sum_lambda_13),
          (((((((sm_dg0_0-sm_dhf0_0)*sm_t)/sm_tref)-sm_dg_0)+sm_dhf0_0)+sm_cpdt_0)-(sm_t*sm_cptdt_0)),
          (((((((sm_dg0_1-sm_dhf0_1)*sm_t)/sm_tref)-sm_dg_1)+sm_dhf0_1)+sm_cpdt_1)-(sm_t*sm_cptdt_1)),
          (((((((sm_dg0_2-sm_dhf0_2)*sm_t)/sm_tref)-sm_dg_2)+sm_dhf0_2)+sm_cpdt_2)-(sm_t*sm_cptdt_2)),
          (((((((sm_dg0_3-sm_dhf0_3)*sm_t)/sm_tref)-sm_dg_3)+sm_dhf0_3)+sm_cpdt_3)-(sm_t*sm_cptdt_3)),
          (((((((sm_dg0_4-sm_dhf0_4)*sm_t)/sm_tref)-sm_dg_4)+sm_dhf0_4)+sm_cpdt_4)-(sm_t*sm_cptdt_4)),
          (((((((sm_dg0_5-sm_dhf0_5)*sm_t)/sm_tref)-sm_dg_5)+sm_dhf0_5)+sm_cpdt_5)-(sm_t*sm_cptdt_5)),
          (((((((sm_dg0_6-sm_dhf0_6)*sm_t)/sm_tref)-sm_dg_6)+sm_dhf0_6)+sm_cpdt_6)-(sm_t*sm_cptdt_6)),
          (((((((sm_dg0_7-sm_dhf0_7)*sm_t)/sm_tref)-sm_dg_7)+sm_dhf0_7)+sm_cpdt_7)-(sm_t*sm_cptdt_7)),
          (((((((sm_dg0_8-sm_dhf0_8)*sm_t)/sm_tref)-sm_dg_8)+sm_dhf0_8)+sm_cpdt_8)-(sm_t*sm_cptdt_8)),
          (((((((sm_dg0_9-sm_dhf0_9)*sm_t)/sm_tref)-sm_dg_9)+sm_dhf0_9)+sm_cpdt_9)-(sm_t*sm_cptdt_9)),
          (((((((sm_dg0_10-sm_dhf0_10)*sm_t)/sm_tref)-sm_dg_10)+sm_dhf0_10)+sm_cpdt_10)-(sm_t*sm_cptdt_10)),
          (((((((sm_dg0_11-sm_dhf0_11)*sm_t)/sm_tref)-sm_dg_11)+sm_dhf0_11)+sm_cpdt_11)-(sm_t*sm_cptdt_11)),
          (((((((sm_dg0_12-sm_dhf0_12)*sm_t)/sm_tref)-sm_dg_12)+sm_dhf0_12)+sm_cpdt_12)-(sm_t*sm_cptdt_12)),
          (((((((sm_dg0_13-sm_dhf0_13)*sm_t)/sm_tref)-sm_dg_13)+sm_dhf0_13)+sm_cpdt_13)-(sm_t*sm_cptdt_13)),
          (((sm_ntx*exp(sm_xl_0))-sm_b_0)+(sm_ntx*exp(sm_xl_1))),
          (((((((((sm_ntx*exp(sm_xl_1))-sm_b_1)+((2*sm_ntx)*exp(sm_xl_2)))+(sm_ntx*exp(sm_xl_4)))+((3*sm_ntx)*exp(sm_xl_5)))+(sm_ntx*exp(sm_xl_7)))+(sm_ntx*exp(sm_xl_9)))+((2*sm_ntx)*exp(sm_xl_11)))+((3*sm_ntx)*exp(sm_xl_13))),
          ((((sm_ntx*exp(sm_xl_3))-sm_b_2)+(sm_ntx*exp(sm_xl_4)))+((2*sm_ntx)*exp(sm_xl_5))),
          (((sm_ntx*exp(sm_xl_6))-sm_b_3/sm_b_3_scale)+(sm_ntx*exp(sm_xl_7))),
          (((sm_ntx*exp(sm_xl_8))-sm_b_4)+(sm_ntx*exp(sm_xl_9))),
          (((sm_ntx*exp(sm_xl_10))-sm_b_5)+(sm_ntx*exp(sm_xl_11))),
          (((sm_ntx*exp(sm_xl_12))-sm_b_6)+((2*sm_ntx)*exp(sm_xl_13))),
          ((sm_ntx*exp(sm_xl_14))-sm_b_7),
          (-0.01+(sm_ntx*exp(sm_xl_15))),
          (((((((((((((((((sm_ntx*exp(sm_xl_0))+(sm_ntx*exp(sm_xl_1)))+(sm_ntx*exp(sm_xl_2)))+(sm_ntx*exp(sm_xl_3)))+(sm_ntx*exp(sm_xl_4)))+(sm_ntx*exp(sm_xl_5)))+(sm_ntx*exp(sm_xl_6)))+(sm_ntx*exp(sm_xl_7)))+(sm_ntx*exp(sm_xl_8)))+(sm_ntx*exp(sm_xl_9)))+(sm_ntx*exp(sm_xl_10)))+(sm_ntx*exp(sm_xl_11)))+(sm_ntx*exp(sm_xl_12)))+(sm_ntx*exp(sm_xl_13)))+(sm_ntx*exp(sm_xl_14)))+(sm_ntx*exp(sm_xl_15)))-sm_ntx),
          (exp(sm_xl_0)-sm_yx_0),
          (exp(sm_xl_1)-sm_yx_1),
          (exp(sm_xl_2)-sm_yx_2),
          (exp(sm_xl_3)-sm_yx_3),
          (exp(sm_xl_4)-sm_yx_4),
          (exp(sm_xl_5)-sm_yx_5),
          (exp(sm_xl_6)-sm_yx_6),
          (exp(sm_xl_7)-sm_yx_7),
          (exp(sm_xl_8)-sm_yx_8),
          (exp(sm_xl_9)-sm_yx_9),
          (exp(sm_xl_10)-sm_yx_10),
          (exp(sm_xl_11)-sm_yx_11),
          (exp(sm_xl_12)-sm_yx_12),
          ((sm_ntx*exp(sm_xl_0))-((sm_y_prime_0*sm_ntx)*(1-exp(sm_xl_3)))),
          ((sm_ntx*exp(sm_xl_1))-((sm_y_prime_1*sm_ntx)*(1-exp(sm_xl_3)))),
          ((sm_ntx*exp(sm_xl_2))-((sm_y_prime_2*sm_ntx)*(1-exp(sm_xl_3)))),
          ((sm_ntx*exp(sm_xl_3))-((sm_y_prime_3*sm_ntx)*(1-exp(sm_xl_3)))),
          ((sm_ntx*exp(sm_xl_4))-((sm_y_prime_4*sm_ntx)*(1-exp(sm_xl_3)))),
          ((sm_ntx*exp(sm_xl_5))-((sm_y_prime_5*sm_ntx)*(1-exp(sm_xl_3)))),
          ((sm_ntx*exp(sm_xl_6))-((sm_y_prime_6*sm_ntx)*(1-exp(sm_xl_3)))),
          ((sm_ntx*exp(sm_xl_7))-((sm_y_prime_7*sm_ntx)*(1-exp(sm_xl_3)))),
          ((sm_ntx*exp(sm_xl_8))-((sm_y_prime_8*sm_ntx)*(1-exp(sm_xl_3)))),
          ((sm_ntx*exp(sm_xl_9))-((sm_y_prime_9*sm_ntx)*(1-exp(sm_xl_3)))),
          ((sm_ntx*exp(sm_xl_10))-((sm_y_prime_10*sm_ntx)*(1-exp(sm_xl_3)))),
          ((sm_ntx*exp(sm_xl_11))-((sm_y_prime_11*sm_ntx)*(1-exp(sm_xl_3)))),
          ((sm_ntx*exp(sm_xl_12))-((sm_y_prime_12*sm_ntx)*(1-exp(sm_xl_3)))),
          ((sm_ntx*exp(sm_xl_13))-((sm_y_prime_13*sm_ntx)*(1-exp(sm_xl_3)))),
          ((sm_ntx*exp(sm_xl_14))-((sm_y_prime_14*sm_ntx)*(1-exp(sm_xl_3)))),
          ((sm_ntx*exp(sm_xl_15))-((sm_y_prime_15*sm_ntx)*(1-exp(sm_xl_3)))),
          ((((0.001*sm_ntx)*exp(sm_xl_0))*sm_mw_0)-sm_m_0),
          ((((0.001*sm_ntx)*exp(sm_xl_1))*sm_mw_1)-sm_m_1),
          ((((0.001*sm_ntx)*exp(sm_xl_2))*sm_mw_2)-sm_m_2),
          ((((0.001*sm_ntx)*exp(sm_xl_3))*sm_mw_3)-sm_m_3),
          ((((0.001*sm_ntx)*exp(sm_xl_4))*sm_mw_4)-sm_m_4),
          ((((0.001*sm_ntx)*exp(sm_xl_5))*sm_mw_5)-sm_m_5),
          ((((0.001*sm_ntx)*exp(sm_xl_6))*sm_mw_6)-sm_m_6),
          ((((0.001*sm_ntx)*exp(sm_xl_7))*sm_mw_7)-sm_m_7),
          ((((0.001*sm_ntx)*exp(sm_xl_8))*sm_mw_8)-sm_m_8),
          ((((0.001*sm_ntx)*exp(sm_xl_9))*sm_mw_9)-sm_m_9),
          ((((0.001*sm_ntx)*exp(sm_xl_10))*sm_mw_10)-sm_m_10),
          ((((0.001*sm_ntx)*exp(sm_xl_11))*sm_mw_11)-sm_m_11),
          ((((0.001*sm_ntx)*exp(sm_xl_12))*sm_mw_12)-sm_m_12),
          ((((0.001*sm_ntx)*exp(sm_xl_13))*sm_mw_13)-sm_m_13),
          ((((0.001*sm_ntx)*exp(sm_xl_14))*sm_mw_14)-sm_m_14),
          ((((0.001*sm_ntx)*exp(sm_xl_15))*sm_mw_15)-sm_m_15),
          (sm_m_0-(sm_ym_0*(((((((((((((((sm_m_0+sm_m_1)+sm_m_2)+sm_m_3)+sm_m_4)+sm_m_5)+sm_m_6)+sm_m_7)+sm_m_8)+sm_m_9)+sm_m_10)+sm_m_11)+sm_m_12)+sm_m_13)+sm_m_14)+sm_m_15))),
          (sm_m_1-(sm_ym_1*(((((((((((((((sm_m_0+sm_m_1)+sm_m_2)+sm_m_3)+sm_m_4)+sm_m_5)+sm_m_6)+sm_m_7)+sm_m_8)+sm_m_9)+sm_m_10)+sm_m_11)+sm_m_12)+sm_m_13)+sm_m_14)+sm_m_15))),
          (sm_m_2-(sm_ym_2*(((((((((((((((sm_m_0+sm_m_1)+sm_m_2)+sm_m_3)+sm_m_4)+sm_m_5)+sm_m_6)+sm_m_7)+sm_m_8)+sm_m_9)+sm_m_10)+sm_m_11)+sm_m_12)+sm_m_13)+sm_m_14)+sm_m_15))),
          (sm_m_3-(sm_ym_3*(((((((((((((((sm_m_0+sm_m_1)+sm_m_2)+sm_m_3)+sm_m_4)+sm_m_5)+sm_m_6)+sm_m_7)+sm_m_8)+sm_m_9)+sm_m_10)+sm_m_11)+sm_m_12)+sm_m_13)+sm_m_14)+sm_m_15))),
          (sm_m_4-(sm_ym_4*(((((((((((((((sm_m_0+sm_m_1)+sm_m_2)+sm_m_3)+sm_m_4)+sm_m_5)+sm_m_6)+sm_m_7)+sm_m_8)+sm_m_9)+sm_m_10)+sm_m_11)+sm_m_12)+sm_m_13)+sm_m_14)+sm_m_15))),
          (sm_m_5-(sm_ym_5*(((((((((((((((sm_m_0+sm_m_1)+sm_m_2)+sm_m_3)+sm_m_4)+sm_m_5)+sm_m_6)+sm_m_7)+sm_m_8)+sm_m_9)+sm_m_10)+sm_m_11)+sm_m_12)+sm_m_13)+sm_m_14)+sm_m_15))),
          (sm_m_6-(sm_ym_6*(((((((((((((((sm_m_0+sm_m_1)+sm_m_2)+sm_m_3)+sm_m_4)+sm_m_5)+sm_m_6)+sm_m_7)+sm_m_8)+sm_m_9)+sm_m_10)+sm_m_11)+sm_m_12)+sm_m_13)+sm_m_14)+sm_m_15))),
          (sm_m_7-(sm_ym_7*(((((((((((((((sm_m_0+sm_m_1)+sm_m_2)+sm_m_3)+sm_m_4)+sm_m_5)+sm_m_6)+sm_m_7)+sm_m_8)+sm_m_9)+sm_m_10)+sm_m_11)+sm_m_12)+sm_m_13)+sm_m_14)+sm_m_15))),
          (sm_m_8-(sm_ym_8*(((((((((((((((sm_m_0+sm_m_1)+sm_m_2)+sm_m_3)+sm_m_4)+sm_m_5)+sm_m_6)+sm_m_7)+sm_m_8)+sm_m_9)+sm_m_10)+sm_m_11)+sm_m_12)+sm_m_13)+sm_m_14)+sm_m_15))),
          (sm_m_9-(sm_ym_9*(((((((((((((((sm_m_0+sm_m_1)+sm_m_2)+sm_m_3)+sm_m_4)+sm_m_5)+sm_m_6)+sm_m_7)+sm_m_8)+sm_m_9)+sm_m_10)+sm_m_11)+sm_m_12)+sm_m_13)+sm_m_14)+sm_m_15))),
          (sm_m_10-(sm_ym_10*(((((((((((((((sm_m_0+sm_m_1)+sm_m_2)+sm_m_3)+sm_m_4)+sm_m_5)+sm_m_6)+sm_m_7)+sm_m_8)+sm_m_9)+sm_m_10)+sm_m_11)+sm_m_12)+sm_m_13)+sm_m_14)+sm_m_15))),
          (sm_m_11-(sm_ym_11*(((((((((((((((sm_m_0+sm_m_1)+sm_m_2)+sm_m_3)+sm_m_4)+sm_m_5)+sm_m_6)+sm_m_7)+sm_m_8)+sm_m_9)+sm_m_10)+sm_m_11)+sm_m_12)+sm_m_13)+sm_m_14)+sm_m_15))),
          (sm_m_12-(sm_ym_12*(((((((((((((((sm_m_0+sm_m_1)+sm_m_2)+sm_m_3)+sm_m_4)+sm_m_5)+sm_m_6)+sm_m_7)+sm_m_8)+sm_m_9)+sm_m_10)+sm_m_11)+sm_m_12)+sm_m_13)+sm_m_14)+sm_m_15))),
          (sm_m_13-(sm_ym_13*(((((((((((((((sm_m_0+sm_m_1)+sm_m_2)+sm_m_3)+sm_m_4)+sm_m_5)+sm_m_6)+sm_m_7)+sm_m_8)+sm_m_9)+sm_m_10)+sm_m_11)+sm_m_12)+sm_m_13)+sm_m_14)+sm_m_15))),
          (sm_m_14-(sm_ym_14*(((((((((((((((sm_m_0+sm_m_1)+sm_m_2)+sm_m_3)+sm_m_4)+sm_m_5)+sm_m_6)+sm_m_7)+sm_m_8)+sm_m_9)+sm_m_10)+sm_m_11)+sm_m_12)+sm_m_13)+sm_m_14)+sm_m_15))),
          (sm_m_15-(sm_ym_15*(((((((((((((((sm_m_0+sm_m_1)+sm_m_2)+sm_m_3)+sm_m_4)+sm_m_5)+sm_m_6)+sm_m_7)+sm_m_8)+sm_m_9)+sm_m_10)+sm_m_11)+sm_m_12)+sm_m_13)+sm_m_14)+sm_m_15))),
          (sm_m_0-((sm_ym_prime_0*(((((((((((((((sm_m_0+sm_m_1)+sm_m_2)+sm_m_3)+sm_m_4)+sm_m_5)+sm_m_6)+sm_m_7)+sm_m_8)+sm_m_9)+sm_m_10)+sm_m_11)+sm_m_12)+sm_m_13)+sm_m_14)+sm_m_15))*(1-sm_ym_3))),
          (sm_m_1-((sm_ym_prime_1*(((((((((((((((sm_m_0+sm_m_1)+sm_m_2)+sm_m_3)+sm_m_4)+sm_m_5)+sm_m_6)+sm_m_7)+sm_m_8)+sm_m_9)+sm_m_10)+sm_m_11)+sm_m_12)+sm_m_13)+sm_m_14)+sm_m_15))*(1-sm_ym_3))),
          (sm_m_2-((sm_ym_prime_2*(((((((((((((((sm_m_0+sm_m_1)+sm_m_2)+sm_m_3)+sm_m_4)+sm_m_5)+sm_m_6)+sm_m_7)+sm_m_8)+sm_m_9)+sm_m_10)+sm_m_11)+sm_m_12)+sm_m_13)+sm_m_14)+sm_m_15))*(1-sm_ym_3))),
          (sm_m_3-((sm_ym_prime_3*(((((((((((((((sm_m_0+sm_m_1)+sm_m_2)+sm_m_3)+sm_m_4)+sm_m_5)+sm_m_6)+sm_m_7)+sm_m_8)+sm_m_9)+sm_m_10)+sm_m_11)+sm_m_12)+sm_m_13)+sm_m_14)+sm_m_15))*(1-sm_ym_3))),
          (sm_m_4-((sm_ym_prime_4*(((((((((((((((sm_m_0+sm_m_1)+sm_m_2)+sm_m_3)+sm_m_4)+sm_m_5)+sm_m_6)+sm_m_7)+sm_m_8)+sm_m_9)+sm_m_10)+sm_m_11)+sm_m_12)+sm_m_13)+sm_m_14)+sm_m_15))*(1-sm_ym_3))),
          (sm_m_5-((sm_ym_prime_5*(((((((((((((((sm_m_0+sm_m_1)+sm_m_2)+sm_m_3)+sm_m_4)+sm_m_5)+sm_m_6)+sm_m_7)+sm_m_8)+sm_m_9)+sm_m_10)+sm_m_11)+sm_m_12)+sm_m_13)+sm_m_14)+sm_m_15))*(1-sm_ym_3))),
          (sm_m_6-((sm_ym_prime_6*(((((((((((((((sm_m_0+sm_m_1)+sm_m_2)+sm_m_3)+sm_m_4)+sm_m_5)+sm_m_6)+sm_m_7)+sm_m_8)+sm_m_9)+sm_m_10)+sm_m_11)+sm_m_12)+sm_m_13)+sm_m_14)+sm_m_15))*(1-sm_ym_3))),
          (sm_m_7-((sm_ym_prime_7*(((((((((((((((sm_m_0+sm_m_1)+sm_m_2)+sm_m_3)+sm_m_4)+sm_m_5)+sm_m_6)+sm_m_7)+sm_m_8)+sm_m_9)+sm_m_10)+sm_m_11)+sm_m_12)+sm_m_13)+sm_m_14)+sm_m_15))*(1-sm_ym_3))),
          (sm_m_8-((sm_ym_prime_8*(((((((((((((((sm_m_0+sm_m_1)+sm_m_2)+sm_m_3)+sm_m_4)+sm_m_5)+sm_m_6)+sm_m_7)+sm_m_8)+sm_m_9)+sm_m_10)+sm_m_11)+sm_m_12)+sm_m_13)+sm_m_14)+sm_m_15))*(1-sm_ym_3))),
          (sm_m_9-((sm_ym_prime_9*(((((((((((((((sm_m_0+sm_m_1)+sm_m_2)+sm_m_3)+sm_m_4)+sm_m_5)+sm_m_6)+sm_m_7)+sm_m_8)+sm_m_9)+sm_m_10)+sm_m_11)+sm_m_12)+sm_m_13)+sm_m_14)+sm_m_15))*(1-sm_ym_3))),
          (sm_m_10-((sm_ym_prime_10*(((((((((((((((sm_m_0+sm_m_1)+sm_m_2)+sm_m_3)+sm_m_4)+sm_m_5)+sm_m_6)+sm_m_7)+sm_m_8)+sm_m_9)+sm_m_10)+sm_m_11)+sm_m_12)+sm_m_13)+sm_m_14)+sm_m_15))*(1-sm_ym_3))),
          (sm_m_11-((sm_ym_prime_11*(((((((((((((((sm_m_0+sm_m_1)+sm_m_2)+sm_m_3)+sm_m_4)+sm_m_5)+sm_m_6)+sm_m_7)+sm_m_8)+sm_m_9)+sm_m_10)+sm_m_11)+sm_m_12)+sm_m_13)+sm_m_14)+sm_m_15))*(1-sm_ym_3))),
          (sm_m_12-((sm_ym_prime_12*(((((((((((((((sm_m_0+sm_m_1)+sm_m_2)+sm_m_3)+sm_m_4)+sm_m_5)+sm_m_6)+sm_m_7)+sm_m_8)+sm_m_9)+sm_m_10)+sm_m_11)+sm_m_12)+sm_m_13)+sm_m_14)+sm_m_15))*(1-sm_ym_3))),
          (sm_m_13-((sm_ym_prime_13*(((((((((((((((sm_m_0+sm_m_1)+sm_m_2)+sm_m_3)+sm_m_4)+sm_m_5)+sm_m_6)+sm_m_7)+sm_m_8)+sm_m_9)+sm_m_10)+sm_m_11)+sm_m_12)+sm_m_13)+sm_m_14)+sm_m_15))*(1-sm_ym_3))),
          (sm_m_14-((sm_ym_prime_14*(((((((((((((((sm_m_0+sm_m_1)+sm_m_2)+sm_m_3)+sm_m_4)+sm_m_5)+sm_m_6)+sm_m_7)+sm_m_8)+sm_m_9)+sm_m_10)+sm_m_11)+sm_m_12)+sm_m_13)+sm_m_14)+sm_m_15))*(1-sm_ym_3))),
          (sm_m_15-((sm_ym_prime_15*(((((((((((((((sm_m_0+sm_m_1)+sm_m_2)+sm_m_3)+sm_m_4)+sm_m_5)+sm_m_6)+sm_m_7)+sm_m_8)+sm_m_9)+sm_m_10)+sm_m_11)+sm_m_12)+sm_m_13)+sm_m_14)+sm_m_15))*(1-sm_ym_3))),
          ((sm_f_add_0-sm_f_in_0)+sm_f_mm_0),
          ((sm_f_add_3-sm_f_in_3)+sm_f_mm_1),
          ((sm_f_add_6-sm_f_in_6)+sm_f_mm_2),
          ((sm_f_add_10-sm_f_in_10)+sm_f_mm_3),
          ((sm_f_add_12-sm_f_in_12)+sm_f_mm_4),
          ((sm_f_add_8-sm_f_in_8)+sm_f_mm_5),
          (sm_f_add_1-sm_f_in_1),
          (sm_f_add_2-sm_f_in_2),
          (sm_f_add_4-sm_f_in_4),
          (sm_f_add_5-sm_f_in_5),
          (sm_f_add_7-sm_f_in_7),
          (sm_f_add_9-sm_f_in_9),
          (sm_f_add_11-sm_f_in_11),
          (sm_f_add_13-sm_f_in_13),
          (sm_f_add_14-sm_f_in_14),
          (sm_f_add_15-sm_f_in_15),
          ((((1000*sm_f_in_kgm_c_tot)*sm_k_c)/12)-sm_f_add_0),
          (((sm_f_in_o2_nm3_lance*sm_p_std)/1000)-(136183*sm_f_add_2)),
          (((1000*((1-sm_cao_fr_dol)*sm_m_dol_in))/sm_mw_9)-sm_f_add_9),
          ((((1000*sm_m_cao_in)/sm_mw_14)+((1000*(sm_cao_fr_dol*sm_m_dol_in))/sm_mw_14))-sm_f_add_14),
          (sm_f_gs_0-sm_f_ex_0),
          (sm_f_gs_1-sm_f_ex_1),
          (sm_f_gs_2-sm_f_ex_2),
          (sm_f_gs_3-sm_f_ex_3),
          (sm_f_gs_4-sm_f_ex_4),
          (sm_f_gs_5-sm_f_ex_5),
          (sm_f_gs_6-sm_f_ex_6),
          (sm_f_gs_7-sm_f_ex_7),
          (sm_f_gs_8-sm_f_ex_8),
          (sm_f_gs_9-sm_f_ex_9),
          (sm_f_gs_10-sm_f_ex_10),
          (sm_f_gs_11-sm_f_ex_11),
          (sm_f_gs_12-sm_f_ex_12),
          (sm_f_gs_13-sm_f_ex_13),
          (sm_f_gs_14-sm_f_ex_14),
          (sm_f_gs_15-sm_f_ex_15),
          (((sm_alpha1*tanh(((sm_beta1*sm_f_in_o2_jetbox1)-sm_teta1)))-sm_o2_switch1)+sm_alpha1),
          (((sm_alpha1*tanh(((sm_beta1*sm_f_in_o2_jetbox2)-sm_teta1)))-sm_o2_switch2)+sm_alpha1),
          (((sm_alpha1*tanh(((sm_beta1*sm_f_in_o2_jetbox3)-sm_teta1)))-sm_o2_switch3)+sm_alpha1),
          ((sm_o2_switch1*sm_f_in_o2_jetbox1)-sm_jetslag1),
          ((sm_o2_switch2*sm_f_in_o2_jetbox2)-sm_jetslag2),
          ((sm_o2_switch3*sm_f_in_o2_jetbox3)-sm_jetslag3),
          (((sm_jetslag1-sm_f_in_o2_nm3)+sm_jetslag2)+sm_jetslag3),
          ((sm_f_in_o2_nm3*sm_bias_o2_sm_star)-sm_f_in_o2_nm3_lance),
          ((sm_ntx*exp(sm_xl_1))-sm_f_gs_1),
          ((sm_k_dc*sm_m_c_float)-sm_f_c_charge),
          (((sm_teta_l*sm_f_in_kgm_c_lance)-sm_f_in_kgm_c_tot)+sm_f_c_charge),
          ((sm_k_cao*sm_m_cao)-sm_m_cao_in),
          ((sm_k_cao*sm_m_dol)-sm_m_dol_in),
          ((sm_cp_0*(sm_t-sm_tref))-(1000*sm_cpdt_0)),
          ((sm_cp_1*(sm_t-sm_tref))-(1000*sm_cpdt_1)),
          ((sm_cp_2*(sm_t-sm_tref))-(1000*sm_cpdt_2)),
          ((sm_cp_3*(sm_t-sm_tref))-(1000*sm_cpdt_3)),
          ((sm_cp_4*(sm_t-sm_tref))-(1000*sm_cpdt_4)),
          ((sm_cp_5*(sm_t-sm_tref))-(1000*sm_cpdt_5)),
          ((sm_cp_6*(sm_t-sm_tref))-(1000*sm_cpdt_6)),
          ((sm_cp_7*(sm_t-sm_tref))-(1000*sm_cpdt_7)),
          ((sm_cp_8*(sm_t-sm_tref))-(1000*sm_cpdt_8)),
          ((sm_cp_9*(sm_t-sm_tref))-(1000*sm_cpdt_9)),
          ((sm_cp_10*(sm_t-sm_tref))-(1000*sm_cpdt_10)),
          ((sm_cp_11*(sm_t-sm_tref))-(1000*sm_cpdt_11)),
          ((sm_cp_12*(sm_t-sm_tref))-(1000*sm_cpdt_12)),
          ((sm_cp_13*(sm_t-sm_tref))-(1000*sm_cpdt_13)),
          ((sm_cp_14*(sm_t-sm_tref))-(1000*sm_cpdt_14)),
          ((sm_cp_15*(sm_t-sm_tref))-(1000*sm_cpdt_15)),
          (((0.0003*sm_t)-(sm_cptdt_13+((7e-008*sm_t)*sm_t)))-0.0454),
          (((0.0003*sm_t)-(sm_cptdt_5+((5e-008*sm_t)*sm_t)))-0.0372),
          (((0.0002*sm_t)-(sm_cptdt_11+((3e-008*sm_t)*sm_t)))-0.0213),
          (((0.0001*sm_t)-(sm_cptdt_7+((2e-008*sm_t)*sm_t)))-0.015),
          (((6e-005*sm_t)-(sm_cptdt_2+((1e-008*sm_t)*sm_t)))-0.009),
          (((4e-005*sm_t)-(sm_cptdt_0+((9e-009*sm_t)*sm_t)))-0.0112),
          (((5e-005*sm_t)-(sm_cptdt_10+((1e-008*sm_t)*sm_t)))-0.0067),
          (((0.0003*sm_t)-(sm_cptdt_9+((1e-007*sm_t)*sm_t)))-0.0627),
          (((0.0002*sm_t)-(sm_cptdt_4+((9e-008*sm_t)*sm_t)))-0.0493),
          (((8e-005*sm_t)-(sm_cptdt_6+((2e-008*sm_t)*sm_t)))-0.0114),
          (((8e-005*sm_t)-(sm_cptdt_3+((2e-008*sm_t)*sm_t)))-0.0114),
          (((0.0001*sm_t)-(sm_cptdt_8+((8e-008*sm_t)*sm_t)))-0.0354),
          (((6e-005*sm_t)-(sm_cptdt_1+((1e-008*sm_t)*sm_t)))-0.0088),
          (((6e-005*sm_t)-(sm_cptdt_12+((1e-008*sm_t)*sm_t)))-0.0078),
          ((((((((((((((((((exp(sm_xl_0)*sm_h_0)+(exp(sm_xl_1)*sm_h_1))+(exp(sm_xl_2)*sm_h_2))+(exp(sm_xl_3)*sm_h_3))+(exp(sm_xl_4)*sm_h_4))+(exp(sm_xl_5)*sm_h_5))+(exp(sm_xl_6)*sm_h_6))+(exp(sm_xl_7)*sm_h_7))+(exp(sm_xl_8)*sm_h_8))+(exp(sm_xl_9)*sm_h_9))+(exp(sm_xl_10)*sm_h_10))+(exp(sm_xl_11)*sm_h_11))+(exp(sm_xl_12)*sm_h_12))+(exp(sm_xl_13)*sm_h_13))+(exp(sm_xl_14)*sm_h_14))+(exp(sm_xl_15)*sm_h_15))*sm_ntx)-sm_e/sm_e_scale),
          (((sm_h_add_0*sm_f_add_0)-sm_h_flow_in_0)+(sm_h_mm_0*sm_f_mm_0)),
          ((sm_h_add_1*sm_f_add_1)-sm_h_flow_in_1),
          ((sm_h_add_2*sm_f_add_2)-sm_h_flow_in_2),
          (((sm_h_add_3*sm_f_add_3)-sm_h_flow_in_3)+(sm_h_mm_1*sm_f_mm_1)),
          ((sm_h_add_4*sm_f_add_4)-sm_h_flow_in_4),
          ((sm_h_add_5*sm_f_add_5)-sm_h_flow_in_5),
          (((sm_h_add_6*sm_f_add_6)-sm_h_flow_in_6)+(sm_h_mm_2*sm_f_mm_2)),
          ((sm_h_add_7*sm_f_add_7)-sm_h_flow_in_7),
          (((sm_h_add_8*sm_f_add_8)-sm_h_flow_in_8)+(sm_h_mm_5*sm_f_mm_5)),
          ((sm_h_add_9*sm_f_add_9)-sm_h_flow_in_9),
          (((sm_h_add_10*sm_f_add_10)-sm_h_flow_in_10)+(sm_h_mm_3*sm_f_mm_3)),
          ((sm_h_add_11*sm_f_add_11)-sm_h_flow_in_11),
          (((sm_h_add_12*sm_f_add_12)-sm_h_flow_in_12)+(sm_h_mm_4*sm_f_mm_4)),
          ((sm_h_add_13*sm_f_add_13)-sm_h_flow_in_13),
          ((sm_h_add_14*sm_f_add_14)-sm_h_flow_in_14),
          ((sm_h_add_15*sm_f_add_15)-sm_h_flow_in_15),
          ((((((((((((((((sm_h_flow_in_0+sm_h_flow_in_1)+sm_h_flow_in_2)+sm_h_flow_in_3)+sm_h_flow_in_4)+sm_h_flow_in_5)+sm_h_flow_in_6)+sm_h_flow_in_7)+sm_h_flow_in_8)+sm_h_flow_in_9)+sm_h_flow_in_10)+sm_h_flow_in_11)+sm_h_flow_in_12)+sm_h_flow_in_13)+sm_h_flow_in_14)+sm_h_flow_in_15)-sm_hflow_in),
          (((((((((((((((((sm_h_0*sm_f_ex_0)+(sm_h_1*sm_f_ex_1))+(sm_h_2*sm_f_ex_2))+(sm_h_3*sm_f_ex_3))+(sm_h_4*sm_f_ex_4))+(sm_h_5*sm_f_ex_5))+(sm_h_6*sm_f_ex_6))+(sm_h_7*sm_f_ex_7))+(sm_h_8*sm_f_ex_8))+(sm_h_9*sm_f_ex_9))+(sm_h_10*sm_f_ex_10))+(sm_h_11*sm_f_ex_11))+(sm_h_12*sm_f_ex_12))+(sm_h_13*sm_f_ex_13))+(sm_h_14*sm_f_ex_14))+(sm_h_15*sm_f_ex_15))-sm_hflow_out),
          ((sm_dhf0_0-sm_h_0)+sm_cpdt_0),
          ((sm_dhf0_1-sm_h_1)+sm_cpdt_1),
          ((sm_dhf0_2-sm_h_2)+sm_cpdt_2),
          ((sm_dhf0_3-sm_h_3)+sm_cpdt_3),
          ((sm_dhf0_4-sm_h_4)+sm_cpdt_4),
          ((sm_dhf0_5-sm_h_5)+sm_cpdt_5),
          ((sm_dhf0_6-sm_h_6)+sm_cpdt_6),
          ((sm_dhf0_7-sm_h_7)+sm_cpdt_7),
          ((sm_dhf0_8-sm_h_8)+sm_cpdt_8),
          ((sm_dhf0_9-sm_h_9)+sm_cpdt_9),
          ((sm_dhf0_10-sm_h_10)+sm_cpdt_10),
          ((sm_dhf0_11-sm_h_11)+sm_cpdt_11),
          ((sm_dhf0_12-sm_h_12)+sm_cpdt_12),
          ((sm_dhf0_13-sm_h_13)+sm_cpdt_13),
          ((sm_dhf0_14-sm_h_14)+sm_cpdt_14),
          ((sm_dhf0_15-sm_h_15)+sm_cpdt_15),
          ((sm_dhf0_0-sm_h_mm_0)+sm_cpdt_mm_0),
          ((sm_dhf0_3-sm_h_mm_1)+sm_cpdt_mm_1),
          ((sm_dhf0_6-sm_h_mm_2)+sm_cpdt_mm_2),
          ((sm_dhf0_10-sm_h_mm_3)+sm_cpdt_mm_3),
          ((sm_dhf0_12-sm_h_mm_4)+sm_cpdt_mm_4),
          ((sm_dhf0_8-sm_h_mm_5)+sm_cpdt_mm_5),
          ((60*sm_q_mm_sm)-sm_q),
          ((sm_cp_0*(sm_t_mm-sm_tref))-(1000*sm_cpdt_mm_0)),
          ((sm_cp_3*(sm_t_mm-sm_tref))-(1000*sm_cpdt_mm_1)),
          ((sm_cp_6*(sm_t_mm-sm_tref))-(1000*sm_cpdt_mm_2)),
          ((sm_cp_10*(sm_t_mm-sm_tref))-(1000*sm_cpdt_mm_3)),
          ((sm_cp_12*(sm_t_mm-sm_tref))-(1000*sm_cpdt_mm_4)),
          ((sm_cp_8*(sm_t_mm-sm_tref))-(1000*sm_cpdt_mm_5)),
          ((((8.314*(1000*sm_f_gs_1))*sm_t_gs)/sm_p_std)-(sm_v_gs*sm_cs_a)),
          ((((sm_sigma_foam/60)*sm_v_gs)*sm_fi_foam)-sm_hf),
          ((1000*((((((((((sm_m_4+sm_m_5)+sm_m_6)+sm_m_7)+sm_m_8)+sm_m_9)+sm_m_10)+sm_m_11)+sm_m_12)+sm_m_13)+sm_m_14))-((sm_hs*sm_cs_a)*sm_dens_slag)),
          (((115*sm_mu)**2)-((sm_sigma_foam**2)*(sm_dens_slag*sm_surf_tension))),
          (((0.75424-sm_surf_tension)-(0.5694*sm_ym_11))-(0.13713*sm_ym_4)),
          (((0.5*tanh(((12.95*sm_hs)-1.289)))+0.5)-sm_fi_foam),
          (((0.5*tanh(((5*sm_hf)-1.24)))+0.5)-sm_e1),
          (((0.5*tanh(((3.24*(1-(sm_m_solid/sm_scrapcharge1st)))-1.292)))+0.5)-sm_e2),
          ((sm_e1*sm_e2)-sm_ef),
          (((1600*exp(((-0.2693*sm_b_mu)-11.6725)))*exp((0.625*sm_b_mu)))-sm_mu),
          ((((((13.8-sm_b_mu)+(39.9355*sm_alpha_mu))-((44.049*sm_alpha_mu)*sm_alpha_mu))+(((30.481-(117.15*sm_alpha_mu))+((129.998*sm_alpha_mu)*sm_alpha_mu))*sm_xg_star))+((((-40.9429+(234.049*sm_alpha_mu))-((300.04*sm_alpha_mu)*sm_alpha_mu))*sm_xg_star)*sm_xg_star))+(((((60.7619-(153.928*sm_alpha_mu))+((211.161*sm_alpha_mu)*sm_alpha_mu))*sm_xg_star)*sm_xg_star)*sm_xg_star)),
          (sm_xm_star-(sm_alpha_mu*(sm_xm_star+sm_xa_star))),
          ((sm_ym_11-sm_xg_star)+sm_ym_p2o5),
          ((((sm_ym_14+sm_ym_9)+sm_ym_4)+sm_ym_7)-sm_xm_star),
          ((sm_ym_13+sm_ym_5)-sm_xa_star),
          ((1000*((((((((((((((((sm_dens_i_0*sm_ym_prime_0)+(sm_dens_i_1*sm_ym_prime_1))+(sm_dens_i_2*sm_ym_prime_2))+(sm_dens_i_3*sm_ym_prime_3))+(sm_dens_i_4*sm_ym_prime_4))+(sm_dens_i_5*sm_ym_prime_5))+(sm_dens_i_6*sm_ym_prime_6))+(sm_dens_i_7*sm_ym_prime_7))+(sm_dens_i_8*sm_ym_prime_8))+(sm_dens_i_9*sm_ym_prime_9))+(sm_dens_i_10*sm_ym_prime_10))+(sm_dens_i_11*sm_ym_prime_11))+(sm_dens_i_12*sm_ym_prime_12))+(sm_dens_i_13*sm_ym_prime_13))+(sm_dens_i_14*sm_ym_prime_14))+(sm_dens_i_15*sm_ym_prime_15)))-sm_dens_slag),
          ((ss_abg*(ss_t-298))-ss_cpdt_1),
          (((60*(((ss_q_power_ss+ss_q_mm_ss)+ss_q_gas)+ss_q_rad_ss))*(ss_t/ss_t_melt))-((ss_n_ex*ss_dh_fus)*ss_k_dm)),
          ((1000*(ss_m_ss/ss_mw))-ss_n),
          ((1000*(ss_m_in_t/ss_mw))-ss_n_in),
          ((((ss_k_cool*ss_n_in)*ss_cpdt_1)/60)-ss_q_chg),
          ((ss_soliddum5*((ss_e_f*ss_q_wall)+(((0.18*ss_p_arc)*ss_k_p))))-ss_q_power_ss),     
          #((ss_soliddum5*((ss_e_f*ss_q_wall)+(1000*((0.18*ss_p_arc)*ss_k_p))))-ss_q_power_ss),
          (((((0.8*(ss_soliddum5*(1-ss_alpha)))*ss_k_p)*ss_p_arc))-ss_q_rad_ss),
          #((1000*(((0.8*(ss_soliddum5*(1-ss_alpha)))*ss_k_p)*ss_p_arc))-ss_q_rad_ss),
          (ss_m_ss-(ss_soliddum5*(ss_m_ss+ss_m_liquid))),
          ((((ss_m_ss*ss_m_liquid)*ss_k_t1)*(ss_t_mm-ss_t))-(ss_q_mm_ss*(ss_m_ss+5))),
          #((1000*(((0.8*ss_alpha)*ss_k_p)*ss_p_arc))-ss_q_wall),
          (((((0.8*ss_alpha)*ss_k_p)*ss_p_arc))-ss_q_wall),
          (1-(((((mm_y_0+mm_y_1)+mm_y_2)+mm_y_3)+mm_y_4)+mm_y_5)),
          (((mm_n_0*mm_mw_0)/1000)-mm_m_mm_0),
          (((mm_n_1*mm_mw_1)/1000)-mm_m_mm_1),
          (((mm_n_2*mm_mw_2)/1000)-mm_m_mm_2),
          (((mm_n_3*mm_mw_3)/1000)-mm_m_mm_3),
          (((mm_n_4*mm_mw_4)/1000)-mm_m_mm_4),
          (((mm_n_5*mm_mw_5)/1000)-mm_m_mm_5),
          ((mm_x_0*mm_f_melt_ss)-mm_melt_0),
          ((mm_x_1*mm_f_melt_ss)-mm_melt_1),
          ((mm_x_2*mm_f_melt_ss)-mm_melt_2),
          ((mm_x_3*mm_f_melt_ss)-mm_melt_3),
          ((mm_x_4*mm_f_melt_ss)-mm_melt_4),
          ((mm_x_5*mm_f_melt_ss)-mm_melt_5),
          (mm_f_sm_0-mm_f_ex_0),
          (mm_f_sm_1-mm_f_ex_1),
          (mm_f_sm_2-mm_f_ex_2),
          (mm_f_sm_3-mm_f_ex_3),
          (mm_f_sm_4-mm_f_ex_4),
          ((((((mm_m_mm_0+mm_m_mm_1)+mm_m_mm_2)+mm_m_mm_3)+mm_m_mm_4)+mm_m_mm_5)-mm_mtot),
          ((((((mm_n_0+mm_n_1)+mm_n_2)+mm_n_3)+mm_n_4)+mm_n_5)-mm_ntot),
          (mm_n_0-(mm_y_0*mm_ntot)),
          (mm_n_1-(mm_y_1*mm_ntot)),
          (mm_n_2-(mm_y_2*mm_ntot)),
          (mm_n_3-(mm_y_3*mm_ntot)),
          (mm_n_4-(mm_y_4*mm_ntot)),
          ((1000*mm_f_in_kgm_c)-(12*mm_f_in)),
          (((mm_beta_i_0*(mm_k_m+mm_k_ml))*(mm_y_0-mm_y_star_c))-mm_f_sm_0),
          (((mm_beta_i_1*mm_k_m)*(mm_y_1-mm_y_slag_1))-mm_f_sm_1),
          (((mm_beta_i_2*mm_k_m)*(mm_y_2-mm_y_slag_2))-mm_f_sm_2),
          (((mm_beta_i_3*mm_k_m)*(mm_y_3-mm_y_slag_3))-mm_f_sm_3),
          (((mm_beta_i_4*mm_k_m)*(mm_y_4-mm_y_slag_4))-mm_f_sm_4),
          (((mm_beta_i_5*mm_k_m)*(mm_y_5-mm_y_slag_5))-mm_f_sm_5),
          ((mm_moltendum5*((mm_e_f*mm_q_wall)+(((0.18*mm_k_p)*mm_p_arc))))-mm_p_power),
          #((mm_moltendum5*((mm_e_f*mm_q_wall)+(1000*((0.18*mm_k_p)*mm_p_arc))))-mm_p_power),
          (((mm_m_slag*mm_k_t2)*(mm_t-mm_t_sm))-mm_q_mm_sm),
          (mm_mtot-(mm_moltendum5*(mm_m_solid+mm_mtot))),
          (((mm_k_mcool*(mm_t-298))/1000)-mm_q_cool),
          (((1-mm_teta_l)*mm_f_in_kgm_c_lance)-mm_f_in_kgm_c),
          ((mm_gamma_d*mm_f_o2_lnc)-(1000*mm_k_ml)),
          (((((0.8*(mm_moltendum5*(1-mm_alpha)))*mm_k_p)*mm_p_arc))-mm_q_rad),
          #((1000*(((0.8*(mm_moltendum5*(1-mm_alpha)))*mm_k_p)*mm_p_arc))-mm_q_rad),
          ((gs_lambda_c+gs_lambda_o)-gs_sum_lambda_0),
          ((gs_lambda_c+(2*gs_lambda_o))-gs_sum_lambda_1),
          ((2*gs_lambda_o)-gs_sum_lambda_2),
          ((gs_lambda_c+(4*gs_lambda_h))-gs_sum_lambda_3),
          ((2*gs_lambda_h)-gs_sum_lambda_4),
          ((gs_lambda_o-gs_sum_lambda_5)+(2*gs_lambda_h)),
          (((9*gs_lambda_c)+(20*gs_lambda_h))-gs_sum_lambda_6),
          ((gs_dg_0+((0.008314*gs_t)*gs_xl_0))+gs_sum_lambda_0),
          ((gs_dg_1+((0.008314*gs_t)*gs_xl_1))+gs_sum_lambda_1),
          ((gs_dg_2+((0.008314*gs_t)*gs_xl_2))+gs_sum_lambda_2),
          ((gs_dg_3+((0.008314*gs_t)*gs_xl_3))+gs_sum_lambda_3),
          ((gs_dg_4+((0.008314*gs_t)*gs_xl_4))+gs_sum_lambda_4),
          ((gs_dg_5+((0.008314*gs_t)*gs_xl_5))+gs_sum_lambda_5),
          ((gs_dg_6+((0.008314*gs_t)*gs_xl_6))+gs_sum_lambda_6),
          (((((((gs_dg0_0-gs_dhf0_0)*gs_t)/gs_tref)-gs_dg_0)+gs_dhf0_0)+gs_cpdt_0)-(gs_tt*gs_cptdt_0)),
          (((((((gs_dg0_1-gs_dhf0_1)*gs_t)/gs_tref)-gs_dg_1)+gs_dhf0_1)+gs_cpdt_1)-(gs_tt*gs_cptdt_1)),
          (((((((gs_dg0_2-gs_dhf0_2)*gs_t)/gs_tref)-gs_dg_2)+gs_dhf0_2)+gs_cpdt_2)-(gs_tt*gs_cptdt_2)),
          (((((((gs_dg0_3-gs_dhf0_3)*gs_t)/gs_tref)-gs_dg_3)+gs_dhf0_3)+gs_cpdt_3)-(gs_tt*gs_cptdt_3)),
          (((((((gs_dg0_4-gs_dhf0_4)*gs_t)/gs_tref)-gs_dg_4)+gs_dhf0_4)+gs_cpdt_4)-(gs_tt*gs_cptdt_4)),
          (((((((gs_dg0_5-gs_dhf0_5)*gs_t)/gs_tref)-gs_dg_5)+gs_dhf0_5)+gs_cpdt_5)-(gs_tt*gs_cptdt_5)),
          (((((((gs_dg0_6-gs_dhf0_6)*gs_t)/gs_tref)-gs_dg_6)+gs_dhf0_6)+gs_cpdt_6)-(gs_tt*gs_cptdt_6)),
          (((((gs_ntx*exp(gs_xl_0))-gs_b_0)+(gs_ntx*exp(gs_xl_1)))+(gs_ntx*exp(gs_xl_3)))+((9*gs_ntx)*exp(gs_xl_6))),
          (((((gs_ntx*exp(gs_xl_0))-gs_b_1)+((2*gs_ntx)*exp(gs_xl_1)))+((2*gs_ntx)*exp(gs_xl_2)))+(gs_ntx*exp(gs_xl_5))),
          ((((((4*gs_ntx)*exp(gs_xl_3))-gs_b_2)+((2*gs_ntx)*exp(gs_xl_4)))+((2*gs_ntx)*exp(gs_xl_5)))+((20*gs_ntx)*exp(gs_xl_6))),
          (((2*gs_ntx)*exp(gs_xl_7))-gs_b_3),
          (((((((((gs_ntx*exp(gs_xl_0))+(gs_ntx*exp(gs_xl_1)))+(gs_ntx*exp(gs_xl_2)))+(gs_ntx*exp(gs_xl_3)))+(gs_ntx*exp(gs_xl_4)))+(gs_ntx*exp(gs_xl_5)))+(gs_ntx*exp(gs_xl_6)))+(gs_ntx*exp(gs_xl_7)))-gs_ntx),
          (exp(gs_xl_0)-gs_yx_0),
          (exp(gs_xl_1)-gs_yx_1),
          (exp(gs_xl_2)-gs_yx_2),
          (exp(gs_xl_3)-gs_yx_3),
          (exp(gs_xl_4)-gs_yx_4),
          (exp(gs_xl_5)-gs_yx_5),
          (exp(gs_xl_6)-gs_yx_6),
          (exp(gs_xl_7)-gs_yx_7),
          ((gs_ntx*exp(gs_xl_0))-((gs_y_dry_0*gs_ntx)*(1-exp(gs_xl_5)))),
          ((gs_ntx*exp(gs_xl_1))-((gs_y_dry_1*gs_ntx)*(1-exp(gs_xl_5)))),
          ((gs_ntx*exp(gs_xl_2))-((gs_y_dry_2*gs_ntx)*(1-exp(gs_xl_5)))),
          ((gs_ntx*exp(gs_xl_3))-((gs_y_dry_3*gs_ntx)*(1-exp(gs_xl_5)))),
          ((gs_ntx*exp(gs_xl_4))-((gs_y_dry_4*gs_ntx)*(1-exp(gs_xl_5)))),
          ((gs_ntx*exp(gs_xl_5))-((gs_y_dry_5*gs_ntx)*(1-exp(gs_xl_5)))),
          ((gs_ntx*exp(gs_xl_6))-((gs_y_dry_6*gs_ntx)*(1-exp(gs_xl_5)))),
          ((gs_ntx*exp(gs_xl_7))-((gs_y_dry_7*gs_ntx)*(1-exp(gs_xl_5)))),
          ((gs_f_in_o2_oxy_burn*gs_bias_o2_gs_star)-gs_f_in_nm3_2),
          (gs_f_in_ch4-gs_f_in_nm3_3),
          (gs_f_h20_electrode-gs_f_in_nm3_5),
          ((1000*gs_f_in_kgm_oil)-gs_f_in_nm3_6),
          (((gs_n_add_0+gs_n_sm_0)+gs_n_ngrs_0)-gs_n_in_0),
          (((gs_n_add_1+gs_n_sm_1)+gs_n_ngrs_1)-gs_n_in_1),
          (((gs_n_add_2+gs_n_sm_2)+gs_n_ngrs_2)-gs_n_in_2),
          (((gs_n_add_3+gs_n_sm_3)+gs_n_ngrs_3)-gs_n_in_3),
          (((gs_n_add_4+gs_n_sm_4)+gs_n_ngrs_4)-gs_n_in_4),
          (((gs_n_add_5+gs_n_sm_5)+gs_n_ngrs_5)-gs_n_in_5),
          (((gs_n_add_6+gs_n_sm_6)+gs_n_ngrs_6)-gs_n_in_6),
          (((gs_n_add_7+gs_n_sm_7)+gs_n_ngrs_7)-gs_n_in_7),
          ((gs_f_in_nm3_0*gs_p_std)-(1000*(136183*gs_n_add_0))),
          ((gs_f_in_nm3_1*gs_p_std)-(1000*(136183*gs_n_add_1))),
          ((gs_f_in_nm3_2*gs_p_std)-(1000*(136183*gs_n_add_2))),
          ((gs_f_in_nm3_3*gs_p_std)-(1000*(136183*gs_n_add_3))),
          ((gs_f_in_nm3_4*gs_p_std)-(1000*(136183*gs_n_add_4))),
          (((gs_f_h2o_star*gs_f_in_nm3_5)*gs_p_std)-(1000*(136183*gs_n_add_5))),
          ((gs_f_in_nm3_7*gs_p_std)-(1000*(136183*gs_n_add_7))),
          (((gs_k_oil*(gs_t_ss/gs_t_melt))*gs_n_oil_gas)-gs_n_add_6),
          ((1000*gs_f_in_kgm_oil)-(128*gs_n_oil_in)),
          (((gs_x_air_0*gs_n_suck)*gs_ea3)-gs_n_ngrs_0),
          (((gs_x_air_1*gs_n_suck)*gs_ea3)-gs_n_ngrs_1),
          (((gs_x_air_2*gs_n_suck)*gs_ea3)-gs_n_ngrs_2),
          (((gs_x_air_3*gs_n_suck)*gs_ea3)-gs_n_ngrs_3),
          (((gs_x_air_4*gs_n_suck)*gs_ea3)-gs_n_ngrs_4),
          (((gs_x_air_5*gs_n_suck)*gs_ea3)-gs_n_ngrs_5),
          (((gs_x_air_6*gs_n_suck)*gs_ea3)-gs_n_ngrs_6),
          (((gs_x_air_7*gs_n_suck)*gs_ea3)-gs_n_ngrs_7),
          (((0.5*gs_n_net)+(0.5*gs_dum7))-gs_n_suck),
          ((((1000*gs_n_net)**2)+1e-006)-((1000*gs_dum7)**2)),
          (((-0.5*gs_n_net)+(0.5*gs_dum8))-gs_n_push),
          ((((1000*gs_n_net)**2)+1e-006)-((1000*gs_dum8)**2)),
          ((gs_n_og-gs_n_net)-((((((((gs_n_add_0+gs_n_sm_0)+(gs_n_add_1+gs_n_sm_1))+(gs_n_add_2+gs_n_sm_2))+(gs_n_add_3+gs_n_sm_3))+(gs_n_add_4+gs_n_sm_4))+(gs_n_add_5+gs_n_sm_5))+(gs_n_add_6+gs_n_sm_6))+(gs_n_add_7+gs_n_sm_7))),
          ((gs_yx_0*(gs_n_og+gs_n_push))-gs_n_ex_0),
          ((gs_yx_1*(gs_n_og+gs_n_push))-gs_n_ex_1),
          ((gs_yx_2*(gs_n_og+gs_n_push))-gs_n_ex_2),
          ((gs_yx_3*(gs_n_og+gs_n_push))-gs_n_ex_3),
          ((gs_yx_4*(gs_n_og+gs_n_push))-gs_n_ex_4),
          ((gs_yx_5*(gs_n_og+gs_n_push))-gs_n_ex_5),
          ((gs_yx_7*(gs_n_og+gs_n_push))-gs_n_ex_7),
          (((1-gs_o2_switch1)*gs_f_in_o2_jetbox1)-gs_jetgas1),
          (((1-gs_o2_switch2)*gs_f_in_o2_jetbox2)-gs_jetgas2),
          (((1-gs_o2_switch3)*gs_f_in_o2_jetbox3)-gs_jetgas3),
          (((gs_jetgas1-gs_f_in_o2_oxy_burn)+gs_jetgas2)+gs_jetgas3),
          (gs_t-(1000*gs_tt)),
          (gs_t_sm-(1000*gs_ttsm)),
          ((gs_cp_0*(gs_tt-gs_ttref))-(1000*gs_cpdt_0)),
          ((gs_cp_1*(gs_tt-gs_ttref))-(1000*gs_cpdt_1)),
          ((gs_cp_2*(gs_tt-gs_ttref))-(1000*gs_cpdt_2)),
          ((gs_cp_3*(gs_tt-gs_ttref))-(1000*gs_cpdt_3)),
          ((gs_cp_4*(gs_tt-gs_ttref))-(1000*gs_cpdt_4)),
          ((gs_cp_5*(gs_tt-gs_ttref))-(1000*gs_cpdt_5)),
          ((gs_cp_6*(gs_tt-gs_ttref))-(1000*gs_cpdt_6)),
          ((gs_cp_7*(gs_tt-gs_ttref))-(1000*gs_cpdt_7)),
          (((4e-005*gs_t)-(gs_cptdt_0+((3e-009*gs_t)*gs_t)))-0.0114),
          (((0.0002*gs_t)-(gs_cptdt_1+((5e-008*gs_t)*gs_t)))-0.0341),
          (((8e-005*gs_t)-(gs_cptdt_2+((2e-008*gs_t)*gs_t)))-0.0172),
          (((0.0002*gs_t)-(gs_cptdt_3+((8e-008*gs_t)*gs_t)))-0.0589),
          (((5e-005*gs_t)-(gs_cptdt_4+((1e-008*gs_t)*gs_t)))-0.0124),
          (((5e-005*gs_t)-(gs_cptdt_5+((4e-009*gs_t)*gs_t)))-0.014),
          (((0.0035*gs_t)-(gs_cptdt_6+((2e-006*gs_t)*gs_t)))-0.8931),
          ((gs_cp_0*(gs_ttsm-gs_ttref))-(1000*gs_cpdt_sm_0)),
          ((gs_cp_1*(gs_ttsm-gs_ttref))-(1000*gs_cpdt_sm_1)),
          ((gs_cp_2*(gs_ttsm-gs_ttref))-(1000*gs_cpdt_sm_2)),
          ((gs_cp_3*(gs_ttsm-gs_ttref))-(1000*gs_cpdt_sm_3)),
          ((gs_cp_4*(gs_ttsm-gs_ttref))-(1000*gs_cpdt_sm_4)),
          ((gs_cp_5*(gs_ttsm-gs_ttref))-(1000*gs_cpdt_sm_5)),
          ((gs_cp_6*(gs_ttsm-gs_ttref))-(1000*gs_cpdt_sm_6)),
          ((gs_cp_7*(gs_ttsm-gs_ttref))-(1000*gs_cpdt_sm_7)),
          (1-(gs_dumgas*gs_tt)),
          (1-(gs_dumgas3*gs_ttsm)),
          ((((((((((gs_h_add_0*gs_n_add_0)+(gs_h_add_1*gs_n_add_1))+(gs_h_add_2*gs_n_add_2))+(gs_h_add_3*gs_n_add_3))+(gs_h_add_4*gs_n_add_4))+(gs_h_add_5*gs_n_add_5))+(gs_h_add_6*gs_n_add_6))+(gs_h_add_7*gs_n_add_7))-gs_fh_in)+((((((((gs_h_sm_0*gs_n_sm_0)+(gs_h_sm_1*gs_n_sm_1))+(gs_h_sm_2*gs_n_sm_2))+(gs_h_sm_3*gs_n_sm_3))+(gs_h_sm_4*gs_n_sm_4))+(gs_h_sm_5*gs_n_sm_5))+(gs_h_sm_6*gs_n_sm_6))+(gs_h_sm_7*gs_n_sm_7))),
          (((((((((gs_h_0*gs_n_ex_0)+(gs_h_1*gs_n_ex_1))+(gs_h_2*gs_n_ex_2))+(gs_h_3*gs_n_ex_3))+(gs_h_4*gs_n_ex_4))+(gs_h_5*gs_n_ex_5))+(gs_h_6*gs_n_ex_6))+(gs_h_7*gs_n_ex_7))-gs_fh_out),
          ((gs_dhf0_0+gs_cpdt_0)-gs_h_0),
          ((gs_dhf0_1+gs_cpdt_1)-gs_h_1),
          ((gs_dhf0_3+gs_cpdt_3)-gs_h_3),
          ((gs_dhf0_4+gs_cpdt_4)-gs_h_4),
          ((gs_dhf0_5+gs_cpdt_5)-gs_h_5),
          ((gs_dhf0_6+gs_cpdt_6)-gs_h_6),
          ((100000*(gs_dhf0_7+gs_cpdt_7))-(100000*gs_h_7)),
          ((100000*(gs_dhf0_2+gs_cpdt_2))-(100000*gs_h_2)),
          ((1000*(gs_dhf0_0+gs_cpdt_sm_0))-(1000*gs_h_sm_0)),
          ((1000*(gs_dhf0_1+gs_cpdt_sm_1))-(1000*gs_h_sm_1)),
          ((1000*(gs_dhf0_2+gs_cpdt_sm_2))-(1000*gs_h_sm_2)),
          ((1000*(gs_dhf0_3+gs_cpdt_sm_3))-(1000*gs_h_sm_3)),
          ((1000*(gs_dhf0_4+gs_cpdt_sm_4))-(1000*gs_h_sm_4)),
          ((1000*(gs_dhf0_5+gs_cpdt_sm_5))-(1000*gs_h_sm_5)),
          ((1000*(gs_dhf0_6+gs_cpdt_sm_6))-(1000*gs_h_sm_6)),
          ((1000*(gs_dhf0_7+gs_cpdt_sm_7))-(1000*gs_h_sm_7)),
          (-(gs_q+(60*(gs_q_cool+gs_q_ss)))),
          (((6.28*rd_r_roof)*(rd_hw_exposedwall+0.001))-rd_a_wall),
          ((rd_m_solid/rd_rho_bulk)-rd_v_scrap),
          ((rd_v_scrap/rd_a_solid)-rd_hs_scrappile),
          ((rd_m_liquid/rd_rho_heel)-rd_v_bath),
          ((rd_v_bath/rd_a_bath)-rd_h_bath),
          (((rd_hs_scrappile-rd_h_furnace_t)+rd_h_bath)+rd_hw_exposedwall),
          ((((rd_a_roof/(rd_a_roof+rd_a_wall))*(1-rd_e_f))*rd_q_wall)-rd_q_radroof),
          ((((rd_a_wall/(rd_a_roof+rd_a_wall))*(1-rd_e_f))*rd_q_wall)-rd_q_radwall),
          ((rd_ua_t1*(rd_t_roof-(0.5*(rd_t_xc_0+rd_t_ec_0))))-rd_q_cw_0),
          ((rd_ua_t2*(rd_t_wall-(0.5*(rd_t_xc_1+rd_t_ec_1))))-rd_q_cw_1),
          (((rd_k_steel*rd_a_wall)/rd_tube)-rd_ua_t2),
          ((((0.001*rd_h_gs)*(rd_a_roof/(rd_a_roof+rd_a_wall)))*(rd_t_gs-rd_t_roof))-rd_q_gascool_0),
          ((((0.001*rd_h_gs)*(rd_a_wall/(rd_a_roof+rd_a_wall)))*(rd_t_gs-rd_t_wall))-rd_q_gascool_1),
          (((sm_ntx*exp(sm_xl_0))+mm_n_0)-(mm_y_c*((sm_ntx*exp(sm_xl_3))+mm_n_1)))
]

for addx in range(len(ode)):
    dae.add_ode(ode[addx])
for addz in range(len(alg)):
    dae.add_alg(alg[addz])
    
#! We will simulate over 60 timesteps.
#dae = {'x':x_all, 'z':z_all, 'p':params, 'ode':ode, 'alg':alg}

#! Initial guess values for algebraic variables (obtained from gproms)
Z_ = [4.039647252239652e+002,
      4.062675089488939e+002,
      4.605567449857415e+000,
      8.806026706100458e+002,
      2.394703645555047e+002,
      2.417731482804334e+002,
      6.009661119046376e+003,
      -1.376143845912923e+002,
      -3.944009972497822e+002,
      -1.299286159999857e-003,
      -5.003734417439769e+001,
      -3.875485800000078e-004,
      -2.281915135921167e+002,
      2.934002543176278e+001,
      -1.075045813533060e+002,
      -4.789572896009366e+000,
      -1.858379156568373e+000,
      -3.352335780496583e+002,
      -9.665510306337200e+001,
      -5.481832490969683e+000,
      -2.437467465921531e+003,
      -1.844027229977790e-001,
      4.016619414990365e+002,
      2.302783724928708e+000,
      1.197351822777523e+002,
      6.012500000000000e-001,
      0.000000000000000e+000,
      8.316008423920995e-003,
      1.559251559253241e-001,
      0.000000000000000e+000,
      0.000000000000000e+000,
      4.161696434961781e-003,
      0.000000000000000e+000,
      8.316008316008318e-001,
      0.000000000000000e+000,
      8.350761759364223e-003,
      1.565767809564281e-001,
      0.000000000000000e+000,
      0.000000000000000e+000,
      4.179088532810238e-003,
      0.000000000000000e+000,
      8.350761651000501e-001,
      4.095638372432791e+002,
      5.100432645000000e+002,
      7.723761926000000e+000,
      1.500000000000000e+002,
      0.000000000000000e+000,
      0.000000000000000e+000,
      3.047293589895977e-001,
      3.794894541817787e-001,
      0.000000000000000e+000,
      1.149348065757172e-002,
      0.000000000000000e+000,
      0.000000000000000e+000,
      1.957471007507650e-008,
      -2.027027027029216e-002,
      0.000000000000000e+000,
      0.000000000000000e+000,
      6.890165062282756e-001,
      0.000000000000000e+000,
      0.000000000000000e+000,
      1.722541265570690e-001,
      0.000000000000000e+000,
      3.445082531141378e+000,
      1.957471007507650e-008,
      0.000000000000000e+000,
      9.734755949475812e-001,
      3.794894541817787e-001,
      0.000000000000000e+000,
      1.837476072146407e-001,
      0.000000000000000e+000,
      3.445082531141378e+000,
      0.000000000000000e+000,
      3.316441955717786e-002,
      6.218328586285410e-001,
      0.000000000000000e+000,
      0.000000000000000e+000,
      1.659693444293164e-002,
      3.316441912681976e+000,
      1.171875000000000e+000,
      7.549516567451065e-014,
      3.312579356866710e+000,
      3.312579356866634e+000,
      3.000000000000000e+002,
      3.312579356866785e+000,
      3.312579356866785e+000,
      5.850911960618273e+002,
      1.691405467907359e-007,
      1.265797883176489e-007,
      9.973790893180734e-008,
      2.070192262000000e+002,
      1.942752246000000e+002,
      1.837968232000000e+002,
      2.070191911846549e+002,
      1.942752000086832e+002,
      1.837968048684892e+002,
      2.980000000000000e-001,
      6.000000000000000e-001,
      1.278380000000468e-004,
      1.249870000000763e-004,
      1.501620000000647e-004,
      -3.514999999953972e-006,
      9.281550000006522e-005,
      1.504599999999634e-004,
      -3.134634999999442e-004,
      9.752900000004061e-005,
      2.535880000000823e-004,
      2.105980000000007e-002,
      4.863920000000022e-003,
      -6.404320000000130e-003,
      1.611960000000079e-003,
      5.447840000000204e-004,
      -2.770800000000007e-002,
      3.355704697986577e+000,
      1.666666666666667e+000,
      6.000000000000000e+002,
      7.849253200000073e-003,
      7.674201800000002e-003,
      9.219946800000045e-003,
      -2.158210000000604e-004,
      5.698871700000052e-003,
      9.238244000000062e-003,
      -1.924665889999999e-002,
      5.988280600000007e-003,
      -3.105898173768614e+001,
      -1.706341199036528e+001,
      -6.387074337251830e-003,
      3.832244602351098e-001,
      0.000000000000000e+000,
      -1.105221507468000e+002,
      -3.935023257982000e+002,
      9.219946800000045e-003,
      -7.452021582099999e+001,
      5.698871700000052e-003,
      -2.418007617560000e+002,
      -2.288792466589000e+002,
      5.988280600000007e-003,
      -1.105298721620000e+002,
      -3.935098750130000e+002,
      1.501620000001758e-004,
      -7.452000351500000e+001,
      9.281550000006522e-005,
      -2.418098495400000e+002,
      -2.288603134635000e+002,
      9.752900000004061e-005,
      1.895195717797149e-002,
      8.642092473154975e+000,
      5.685587153391403e-003,
      4.586373637069087e-001,
      3.506112077924715e-001,
      0.000000000000000e+000,
      0.000000000000000e+000,
      7.430181069548107e-001,
      1.842383528236292e+000,
      5.780874512664180e-004,
      1.256399891908178e-003,
      5.059172224208453e-002,
      5.404500000000034e-003,
      9.103550000000000e+000,
      8.100298659999949e-003,
      9.001159599999475e-004,
      5.802150525740002e-002,
      2.442490654175344e-015,
      7.430181069548107e-001,
      1.842383528236292e+000,
      5.780874512664180e-004,
      1.256399891908178e-003,
      5.059172224208453e-002,
      3.441691376337985e-015,
      2.714440487399439e-003,
      9.832306654358262e-001,
      8.893653133814627e-004,
      1.932922910627966e-004,
      1.297223647232937e-002,
      8.881784197001252e-016,
      1.657800206300001e+002,
      0.000000000000000e+000,
      2.256643319997332e-001,
      9.066000000000019e-003,
      6.417512004372963e+000,
      0.000000000000000e+000,
      2.182622242493801e-001,
      2.832486895214870e-003,
      2.746216161141752e-001,
      3.740785459172002e-012,
      0.000000000000000e+000,
      0.000000000000000e+000,
      0.000000000000000e+000,
      0.000000000000000e+000,
      6.000000000000000e+002,
      0.000000000000000e+000,
      9.475978588985717e+000,
      5.398240000000000e+001,
      1.452851852757493e-001,
      9.175976419877403e+000,
      5.422598344535068e-003,
      1.039175634431899e-009,
      1.298969545593387e-006,
      0.000000000000000e+000,
      0.000000000000000e+000,
      3.026792270090466e-001,
      6.102737745973801e+000,
      5.998044444444444e+001,
      5.398240000000000e+001,
      1.865434801840057e+000,
      1.025248761997475e+000,
      3.188597115089675e-002,
      0.000000000000000e+000,
      0.000000000000000e+000,
      -1.286144000000000e+001,
      -6.102737745973801e-001,
      3.051368872986900e-001,
      9.175976419877403e+000,
      2.980000000000000e+002,
      5.422598344535068e-003,
      0.000000000000000e+000,
      0.000000000000000e+000,
      0.000000000000000e+000,
      3.496717839603613e+001,
      2.656685207755980e+002,
      4.614026847591236e+002,
      9.958673962561486e+000,
      2.406600163421233e+002,
      7.120213750638085e+002,
      1.347653247373219e+002,
      3.654666671168837e+002,
      3.482271070215502e+002,
      5.789284494011121e+002,
      3.195677740910293e+002,
      7.809704588501529e+002,
      4.096177928999046e+002,
      1.511339612938495e+003,
      3.496717839603613e+001,
      2.307013423795618e+002,
      9.958673962561486e+000,
      3.482271070215502e+002,
      1.347653247373219e+002,
      3.195677740910293e+002,
      4.096177928999046e+002,
      -5.702230000000075e+000,
      -1.708367648122867e+002,
      -4.842279999999985e+000,
      -3.511859999999988e+000,
      -2.336231423549488e+002,
      -6.934049808873720e+002,
      -3.511859999999988e+000,
      -3.442302485665529e+002,
      1.305316999999999e+001,
      -5.666300437883959e+002,
      -3.469599999999964e+000,
      -7.666494310921501e+002,
      -5.009679999999973e+000,
      -1.489577613993174e+003,
      -5.866600191651844e+000,
      -1.901045544930464e+001,
      -9.152441760065825e+001,
      -1.292361070195152e+000,
      -1.410647499634047e+000,
      -3.731936928962478e+000,
      -2.631173617539130e+001,
      -4.257160322013252e+000,
      -7.242407926821228e+001,
      -2.465400852521063e+000,
      -6.336664543561648e+001,
      -2.870865960629228e+000,
      -8.110979730973953e+001,
      -4.362520837406945e+000,
      -1.261428048195127e+000,
      -5.866598234183218e+000,
      6.000000000000000e+002,
      7.430181069548107e-001,
      0.000000000000000e+000,
      5.798863889072692e-008,
      1.842383528236292e+000,
      0.000000000000000e+000,
      0.000000000000000e+000,
      5.780874512664180e-004,
      0.000000000000000e+000,
      3.441691376337985e-015,
      0.000000000000000e+000,
      1.256399891908178e-003,
      0.000000000000000e+000,
      5.059172224208453e-002,
      0.000000000000000e+000,
      0.000000000000000e+000,
      0.000000000000000e+000,
      0.000000000000000e+000,
      5.798863889072692e-008,
      0.000000000000000e+000,
      0.000000000000000e+000,
      0.000000000000000e+000,
      1.957471007507650e-008,
      -2.027027027029216e-002,
      0.000000000000000e+000,
      0.000000000000000e+000,
      0.000000000000000e+000,
      0.000000000000000e+000,
      0.000000000000000e+000,
      0.000000000000000e+000,
      0.000000000000000e+000,
      0.000000000000000e+000,
      0.000000000000000e+000,
      0.000000000000000e+000,
      0.000000000000000e+000,
      0.000000000000000e+000,
      0.000000000000000e+000,
      1.957471007507650e-008,
      -2.027027027029216e-002,
      7.430181069548107e-001,
      1.842383528236292e+000,
      5.780874512664180e-004,
      1.256399891908178e-003,
      5.059172224208453e-002,
      3.441691376337985e-015,
      3.530459544647022e+000,
      2.832486895214870e-003,
      5.544521841649441e-009,
      0.000000000000000e+000,
      2.746216161141752e-001,
      2.439852516351366e-001,
      2.394640839592754e-002,
      3.740785459172002e-012,
      1.416246219488204e-002,
      0.000000000000000e+000,
      8.497477319173674e-002,
      0.000000000000000e+000,
      5.664984879449109e-002,
      0.000000000000000e+000,
      3.904840505504770e-003,
      7.643627175646373e-009,
      0.000000000000000e+000,
      3.785908461223196e-001,
      3.363558344930501e-001,
      3.301229941213246e-002,
      5.157096971686315e-012,
      1.952424074041781e-002,
      0.000000000000000e+000,
      1.171454444734486e-001,
      0.000000000000000e+000,
      7.809696298229907e-002,
      0.000000000000000e+000,
      1.757181667101715e-002,
      3.904848149114950e-001,
      3.904848149114848e-003,
      1.200997645355306e-004,
      5.482876375140222e-010,
      0.000000000000000e+000,
      5.414883724762754e-002,
      6.189015733890457e-002,
      1.350132962520245e-002,
      7.255307465925398e-013,
      3.547499988053704e-003,
      0.000000000000000e+000,
      1.208999996248006e-002,
      0.000000000000000e+000,
      1.201599996270963e-002,
      0.000000000000000e+000,
      4.588199985760988e-003,
      5.607999982596179e-002,
      2.800999991307007e-004,
      5.502544700465295e-004,
      2.512059227299801e-009,
      0.000000000000000e+000,
      2.480907423804067e-001,
      2.835587218619685e-001,
      6.185829761258288e-002,
      3.324118758030181e-012,
      1.625338512082797e-002,
      0.000000000000000e+000,
      5.539208630379566e-002,
      0.000000000000000e+000,
      5.505304458448335e-002,
      0.000000000000000e+000,
      2.102150292630867e-002,
      2.569386435001524e-001,
      1.283318724043858e-003,
      0.000000000000000e+000,
      0.000000000000000e+000,
      7.793817273515913e-005,
      0.000000000000000e+000,
      1.691405467907359e-007,
      1.265797883176489e-007,
      9.973790893180734e-008,
      3.501534511562987e-005,
      2.459131680521143e-005,
      1.833151081431783e-005,
      7.793817273515913e-005,
      0.000000000000000e+000,
      3.376999999999997e-002,
      8.721870000000001e+000,
      9.197720000000000e+000,
      1.412814000000000e+001,
      1.497853000000000e+001,
      4.623419999999999e+001,
      1.412814000000000e+001,
      1.862569000000000e+001,
      1.053317000000000e+001,
      1.862569000000000e+001,
      8.350400000000001e+000,
      2.646033000000000e+001,
      9.750320000000000e+000,
      5.651870000000000e+001,
      1.547894000000000e+001,
      9.338940000000001e+000,
      9.560000000000124e-003,
      2.359999999999995e-002,
      2.339999999999998e-002,
      2.939999999999998e-002,
      3.830000000000000e-002,
      1.247999999999999e-001,
      2.939999999999998e-002,
      3.779999999999995e-002,
      -4.199999999999982e-003,
      8.129999999999993e-002,
      1.969999999999994e-002,
      8.790000000000009e-002,
      2.459999999999996e-002,
      1.093999999999998e-001,
      1.809000000000000e+003,
      1.488591050424759e-008,
      2.980000000000000e+002,
      1.816579109359395e-010,
      1.002707865215192e+001,
      1.422229171547773e-003,
      3.580011940632761e+003,
      7.302232509768258e-002,
      6.840083888846635e-001,
      4.314689350698003e+000,
      7.727220226618425e-002,
      7.017527889079067e-002,
      5.422598344535068e-003,
      5.398240000000000e+001,
      1.618433067424840e+001,
      8.807523725301911e-001,
      5.605304458448335e-002,
      6.121428367867445e-001,
      8.287980053889155e-002,
      7.318096757960024e-004,
      3.340907439586260e-009,
      0.000000000000000e+000,
      3.299477162521134e-001,
      3.771182745636932e-001,
      8.226830164109811e-002,
      4.420908084057373e-012,
      2.161615242281134e-002,
      0.000000000000000e+000,
      7.366857867817300e-002,
      0.000000000000000e+000,
      7.321767091785947e-002,
      0.000000000000000e+000,
      2.795749980903139e-002,
      3.417149621399436e-001,
      1.706746806265813e-003,
      1.353985991998400e+001,
      2.256643319997332e-001,
      3.376999999999997e-002,
      -1.018781300000000e+002,
      9.197720000000000e+000,
      1.412814000000000e+001,
      -2.555714700000000e+002,
      -7.848458000000001e+002,
      1.412814000000000e+001,
      -3.667243100000000e+002,
      1.053317000000000e+001,
      -5.836043100000001e+002,
      8.350400000000001e+000,
      -8.249196700000000e+002,
      9.750320000000000e+000,
      -1.614381300000000e+003,
      -6.194210600000000e+002,
      9.338940000000001e+000,
      1.667600000000000e-001,
      6.976632000000001e+001,
      6.976632000000001e+001,
      4.123520000000000e+001,
      4.814816000000000e+001,
      5.201396000000000e+001,
      1.311882617739881e+002,
      -1.864422645053294e-001,
      1.239056995157842e-001,
      0.000000000000000e+000,
      8.686698138582472e-009,
      1.285363187936622e+002,
      0.000000000000000e+000,
      0.000000000000000e+000,
      4.033103411303729e-002,
      0.000000000000000e+000,
      1.789679515695752e-013,
      0.000000000000000e+000,
      5.180790082281206e-002,
      0.000000000000000e+000,
      2.435898337187445e+000,
      0.000000000000000e+000,
      0.000000000000000e+000,
      0.000000000000000e+000,
      1.667600000000000e-001,
      6.976632000000001e+001,
      6.976632000000001e+001,
      4.123520000000000e+001,
      4.814816000000000e+001,
      5.201396000000000e+001,
      4.798000000000002e-002,
      9.475978588985717e+000,
      9.175976419877403e+000,
      9.674265232974910e+002,
      0.000000000000000e+000,
      0.000000000000000e+000,
      -6.387074337251830e-003,
      6.417512004372963e+000,
      0.000000000000000e+000,
      0.000000000000000e+000,
      1.809000000000000e+003,
      8.547148147242507e-001,
      2.980000000000000e+002,
      5.422598344535068e-003,
      0.000000000000000e+000,
      0.000000000000000e+000,
      0.000000000000000e+000]#mm_y_c

# Row vector of parameter values 
# Had to divide the contents into 2 parts because python not taking more than
# 255 arguments for a function         
par_val1 = vertcat(-1.292,
10,
32.1536,
0.9,
1,##
3.14,
1809,
293,
0.008314,
53.9824,
3.98802,
0.37,
32.1536,
0.0961354,
0.6,
5,
0.693895,
0.535233,
0.45,
0.01,
11.6865,
3.2,
0.0113723,
3.14,
74.986,
36.91,
1809,
32.1536,
900,
-1.289,
0.9706,
2.2,
7.82,
0.02399,
8.314,
55.8,
-1.24,
0.293,
0.8,
0.535233,
70.7379,
0.001,
53.9824,
0.002,
0.7,
293,
2,
4.16759,
3.41297,
293,
3.24,
0.96,
8.95,
0.293,
298,
0.78*0.509895,
0.006,
0.75,
0.039,
0.65,
13,##
1.3,
1.2,
0.134,
0.436634,
32.1536,
0.5822,
9,
101325,
0.8,
0.000855181,
1.60768,
14.96,
0.827393,
130,
101325,
0.7,
6.4,
0.000506402,
0.0000864619,
0.15,
0.75,
0.0005,
0.78*0.509895,#
3.14,
0.985,
12.95,
0.11,
28.41,
29.96,
46.02,
48.79,
150.6,
46.02,
60.67,
34.31,
60.67,
27.2,
86.19,
31.76,
184.1,
50.42,
30.42,
300,
298,
40.3,
12.01,
28.09,
60.08,
26.98,
101.96,
56.08,
28.01,
32,
55.85,
71.85,
159.7,
54.94,
70.95,
24.3,
28.01,
0.002,
0.912,
0.0006,
0.0484,
0.037,
0,
1000,
1000,
1000,
1000,
15.8,
38.4,
1000,
15.6,
1000,
16.1,
1000,
20.35,
1000,
29,
20.7,
1000,
0.00055,
0.14205,
0.1498,
0.2301,
0.24395,
0.753,
0.2301,
0.30335,
0.17155,
0.30335,
0.136,
0.43095,
0.1588,
0.9205,
0.2521,
0.1521,
28,
44,
32,
16,
2,
18,
128.26,
28,
0,
0,
0.16,
0,
0,
0.04,
0,
0.8,
-110.53,
-393.51,
0,
-74.52,
0)

mis_k_p = 0.7 #1.1

# modified params part
par_val2 = vertcat(-1.292,
10,
32.1536,
0.9,
1.00, #1.05,## Original: 1.00
3.14,
1809,
293,
0.008314,
53.9824,
3.98802,
0.37,
32.1536,
0.0961354,
0.6,
5,
0.693895,
0.535233,
0.45,
0.01,
11.6865,
3.2,
0.0113723,
3.14,
74.986,
36.91,
1809,
32.1536,
900,
-1.289,
0.9706,
2.2,
7.82,
0.02399,
8.314,
55.8,
-1.24,
0.293,
0.8,
0.535233,
70.7379,
0.001,
53.9824,
0.002,
0.7,
293,
2,
4.16759,
3.41297,
293,
3.24,
0.96,
8.95,
0.293,
298,
mis_k_p*0.509895, #Original: 0.509895
0.006,
0.75,
0.039,
0.65,
12.35, #13,##Original : 12.35
1.3,
1.2,
0.134,
0.436634,
32.1536,
0.5822,
9,
101325,
0.8,
0.000855181,
1.60768,
14.96,
0.827393,
130,
101325,
0.7,
6.4,
0.000506402,
0.0000864619,
0.15,
0.75,
0.0005,
mis_k_p*0.509895, # 0.53538975, # Original param: 0.509895
3.14,
0.985,
12.95,
0.11,
28.41,
29.96,
46.02,
48.79,
150.6,
46.02,
60.67,
34.31,
60.67,
27.2,
86.19,
31.76,
184.1,
50.42,
30.42,
300,
298,
40.3,
12.01,
28.09,
60.08,
26.98,
101.96,
56.08,
28.01,
32,
55.85,
71.85,
159.7,
54.94,
70.95,
24.3,
28.01,
0.002,
0.912,
0.0006,
0.0484,
0.037,
0,
1000,
1000,
1000,
1000,
15.8,
38.4,
1000,
15.6,
1000,
16.1,
1000,
20.35,
1000,
29,
20.7,
1000,
0.00055,
0.14205,
0.1498,
0.2301,
0.24395,
0.753,
0.2301,
0.30335,
0.17155,
0.30335,
0.136,
0.43095,
0.1588,
0.9205,
0.2521,
0.1521,
28,
44,
32,
16,
2,
18,
128.26,
28,
0,
0,
0.16,
0,
0,
0.04,
0,
0.8,
-110.53,
-393.51,
0,
-74.52,
0)

par_val = vertcat(par_val1,-241.81,
-228.86,
0,
80,
0.2,
0.05,
0.5,
0.3,
0.3,
25.5676,
24.9974,
30.0324,
-0.703,
18.5631,
30.092,
-62.6927,
19.5058,
0.11,
46.02,
34.31,
27.2,
31.76,
46.02,
0.00011,
0.04602,
0.0272,
0.03176,
0.04602,
0.03431,
312,
302,
0,
-137.36,
0,
0,
-248.61,
-749.86,
0,
-363.29,
0,
-570.12,
0,
-797.17,
0,
-1577.9,
0,
0,
12.01,
55.85,
54.94,
28.09,
26.98,
24.3,
0.00011,
0.04602,
0.0272,
0.03176,
0.04602,
0.03431,
0,
-110.6,
0,
0,
-270.55,
-831.08,
0,
-385.35,
0,
-602.23,
0,
-851.38,
0,
-1670.9,
-634.9,
0,
-110.53,
-393.51,
0,
-74.52,
0,
-241.81,
-228.86,
0,
0,
0,
0,
0,
4.54747,
4.15885,
0.05494,
4.54808,
0.0243,
2.50311,
0.02809,
2.95233,
0.02698,
3.51586,
2.70918,
0,
-137.16,
-394.38,
0,
-50.45,
0,
-228.42,
25,
0,
0.00055,
-110.458,
0.1498,
0.2301,
-270.306,
-830.327,
0.2301,
-385.047,
0.17155,
-601.927,
0.136,
-850.949,
0.1588,
-1669.98,
-634.648,
0.1521,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0.01,
0.0000000000001,
0)

par_val_mismatch = vertcat(par_val2,-241.81,
-228.86,
0,
80,
0.2,
0.05,
0.5,
0.3,
0.3,
25.5676,
24.9974,
30.0324,
-0.703,
18.5631,
30.092,
-62.6927,
19.5058,
0.11,
46.02,
34.31,
27.2,
31.76,
46.02,
0.00011,
0.04602,
0.0272,
0.03176,
0.04602,
0.03431,
312,
302,
0,
-137.36,
0,
0,
-248.61,
-749.86,
0,
-363.29,
0,
-570.12,
0,
-797.17,
0,
-1577.9,
0,
0,
12.01,
55.85,
54.94,
28.09,
26.98,
24.3,
0.00011,
0.04602,
0.0272,
0.03176,
0.04602,
0.03431,
0,
-110.6,
0,
0,
-270.55,
-831.08,
0,
-385.35,
0,
-602.23,
0,
-851.38,
0,
-1670.9,
-634.9,
0,
-110.53,
-393.51,
0,
-74.52,
0,
-241.81,
-228.86,
0,
0,
0,
0,
0,
4.54747,
4.15885,
0.05494,
4.54808,
0.0243,
2.50311,
0.02809,
2.95233,
0.02698,
3.51586,
2.70918,
0,
-137.16,
-394.38,
0,
-50.45,
0,
-228.42,
25,
0,
0.00055,
-110.458,
0.1498,
0.2301,
-270.306,
-830.327,
0.2301,
-385.047,
0.17155,
-601.927,
0.136,
-850.949,
0.1588,
-1669.98,
-634.648,
0.1521,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0.01,
0.0000000000001,
0)

#%%##########################################
# State variable declarations
#############################################

dae.set_start('sm_b_0', 0.01)
dae.set_start('sm_b_1', 2)
dae.set_start('sm_b_2', 2)
#dae.set_start('sm_b_3', 0.05)
dae.set_start('sm_b_3', 0.05*sm_b_3_scale)
dae.set_start('sm_b_4', 0.3)
dae.set_start('sm_b_5', 0.2)
dae.set_start('sm_b_6', 0.09)
dae.set_start('sm_b_7', 1)
dae.set_start('sm_m_cao', 0)
dae.set_start('sm_m_c_float', 0)
dae.set_start('sm_m_dol', 0)
dae.set_start('sm_e', -1323.1746*sm_e_scale)
dae.set_start('ss_m_ss', 53.9824)
dae.set_start('ss_t', 300)
dae.set_start('mm_n_0', 0.45)
dae.set_start('mm_n_1', 163)
dae.set_start('mm_n_2', 0.147439)
dae.set_start('mm_n_3', 0.032044)
dae.set_start('mm_n_4', 2.15053763)
dae.set_start('mm_t', 1809)
dae.set_start('gs_b_0', 0.005)
dae.set_start('gs_b_1', 0.2)
dae.set_start('gs_b_2', 0.005)
dae.set_start('gs_b_3', 1)
dae.set_start('gs_n_oil_gas', 0)
dae.set_start('gs_t', 298)
dae.set_start('rd_t_roof', 298)
dae.set_start('rd_t_wall', 298)

# Disturbance state addition
if disturb_add ==1:
    dae.set_start('d_sm_b_3',0)   

for i4_ in range(len(dae.z)):
    dae.set_guess(dae.z[i4_],[Z_[i4_]])  

#%%##########################################
# Initial values for state variables
#############################################
       
# Initial value for states
if disturb_add == 1:       
    X0 = [0.01, 2, 2, 0.05*sm_b_3_scale, 0.3, 0.2, 0.09, 1,  # sm_b          (Moles of C, O, Fe, Mn, Mg, Si, Al, CaO)
          0,                                    # sm_m_cao      (Amount of lime floating on top of slag)
          0,                                    # sm_m_c_float  (Amount of carbon floating on top of slag)
          0,                                    # sm_m_dol      (Amount of Dolomite floating on top of slag)
          -1323.1746*sm_e_scale,                           # sm_e          (Energy holdup of slag)
          53.9824,                              # ss_m_ss       (Solid scrap charged)
          300,                                  # ss_t          (Scrap temperature)
          0.45,                                 # mm_n[1]       (Moles of C)
          163,                                  # mm_n[2]       (Moles of Fe)
          0.147439,                             # mm_n[3]       (Moles of Mn)
          0.032044,                             # mm_n[4]       (Moles of Si)
          2.15053763,                           # mm_n[5]       (Moles of Al)
          1809,                                 # mm_t          (Molten metal temperature)
          0.005,                                # gs_b[1]       (Moles of C)
          0.2,                                  # gs_b[2]       (Moles of O)
          0.005,                                # gs_b[3]       (Moles of H)
          1,                                    # gs_b[4]       (Moles of N)
          0,                                    # gs_n_oil_gas  (Moles of volatile material)
          298,                                  # gs_t          (Gas zone temperature)
          298,                                  # rd_t_roof     (Furnace roof temperature)
          298,                                  # rd_t_wall     (Furnace wall temperature)
          0]                                    # d_sm_b_3      (Disturbance state for sm_b_3)
else:
    X0 = [0.01, 2, 2, 0.05*sm_b_3_scale, 0.3, 0.2, 0.09, 1,  # sm_b          (Moles of C, O, Fe, Mn, Mg, Si, Al, CaO)
          0,                                    # sm_m_cao      (Amount of lime floating on top of slag)
          0,                                    # sm_m_c_float  (Amount of carbon floating on top of slag)
          0,                                    # sm_m_dol      (Amount of Dolomite floating on top of slag)
          -1323.1746*sm_e_scale,                           # sm_e          (Energy holdup of slag)
          53.9824,                              # ss_m_ss       (Solid scrap charged)
          300,                                  # ss_t          (Scrap temperature)
          0.45,                                 # mm_n[1]       (Moles of C)
          163,                                  # mm_n[2]       (Moles of Fe)
          0.147439,                             # mm_n[3]       (Moles of Mn)
          0.032044,                             # mm_n[4]       (Moles of Si)
          2.15053763,                           # mm_n[5]       (Moles of Al)
          1809,                                 # mm_t          (Molten metal temperature)
          0.005,                                # gs_b[1]       (Moles of C)
          0.2,                                  # gs_b[2]       (Moles of O)
          0.005,                                # gs_b[3]       (Moles of H)
          1,                                    # gs_b[4]       (Moles of N)
          0,                                    # gs_n_oil_gas  (Moles of volatile material)
          298,                                  # gs_t          (Gas zone temperature)
          298,                                  # rd_t_roof     (Furnace roof temperature)
          298]#,                                  # rd_t_wall     (Furnace wall temperature)
          #-0.001]                                    # d_sm_b_3      (Disturbance state for sm_b_3)

#%%##########################################
# Matrix of inputs (Total 12)
#############################################

# Time:  0          1           2           3               4               5               6               7           8           9           10          11          12          13          14          15          16              17          18          19          20          21          22              23          24          25          26          27          28          29              30          31          32          33          34              35          36          37          38              39              40              41              42              43              44              45              46              47              48              49              50              51              52              53              54              55              56              57              58              59
# gs_f_in_kgm_oil       (volatile input)
U_dm = DM([[0.15,      0,          0,          0,              0,              0,              0,              0,          0,          0,          0,          0,          0,          0,          0,          0,          0,              0,          0,          0,          0,          0,          0,              0,          0,          0.15,       0,          0,          0,          0,              0,          0,          0,          0,          0,              0,          0,          0,          0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0],

# sm_f_in_cao           (Flow rate of lime)
        [0.7983226, 0,          0,          0,              0,              0,              0,              0,          0.01238307, 0,          0,          0,          0,          0,          0,          0,          0,              0.27210885, 0.14912872, 0.16342874, 0.14708586, 0.08375723, 0,              0,          0,          0,          0,          0,          0,          0,              0,          0,          0,          0,          0,              0,          0,          0,          0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0],

# sm_f_in_kgm_c_charge  (Carbon charged additions)
        [1.102229,  0,          0,          0,              0,              0,              0,              0,          0,          0,          0,          0,          0,          0,          0,          0,          0,              0,          0,          0,          0,          5.00E-04,   0.029,          0.04,       0.0075,     0,          0,          0,          0,          0,              0,          0,          0,          0,          0,              0,          0,          0,          0.00472973,     0.017792793,    0.018693693,    0.02027027,     0.022972973,    0.023198199,    0.023873873,    0.012387387,    0.011711712,    0.011486487,    0.011711712,    0.013738738,    0.010135135,    0.011036036,    0.012162162,    0.012612613,    0.015990991,    0.016441442,    0.016666668,    0.017117117,    0.017117117,    0.002927928],

# sm_f_in_kgm_c_lance   (Carbon lanced additions)
        [0,         0,          0,          0,              0,              0,              0,              0,          0,          0,          0,          0,          0,          0,          0,          0,          0,              0,          0,          0,          0,          9.53E-03,   0.026535425,    0.03538057, 0,          0,          0,          0,          0,          0,              0,          0,          0,          0,          0,              0,          0,          0,          0.018370679,    0.037194956,    0.037421755,    0.04105053,     0.0460401,      0.049668875,    0.05012247,     0.023587046,    0.022453053,    0.023133447,    0.023587046,    0.028349815,    0.020638665,    0.020411866,    0.025401434,    0.025401434,    0.032658987,    0.03311258,     0.034246575,    0.03333938,     0.03515377,     0],

# gs_f_h20_electrode    (Water from sprays)
        [7.723762,  435.02954,  408.90503,  408.22354,      408.33713,      397.66016,      397.09222,      403.79373,  0,          0,          0,          0,          0,          0,          0,          0,          298.84143,      423.78464,  415.0386,   412.99408,  4.09E+02,   406.9741,   396.29715,      409.70013,  492.617,    580.3044,   578.9414,   578.9414,   393.79828,  417.76465,      418.9005,   417.4239,   411.06314,  402.99863,  394.82053,      404.47525,  415.49295,  416.7424,   417.76465,      408.90503,      395.50204,      398.56882,      408.10995,      417.65106,      419.12766,      411.85825,      397.88733,      396.29715,      409.1322,       417.87823,      420.26352,      414.4707,       402.08997,      395.1613,       402.5443,       412.53976,      419.24124,      419.6956,       415.60654,      1.9309405],

# ss_m_in_t             (2nd charge addition)
        [0,         0,          0,          0,              0,              0,              0,              0,          0,          0,          0,          0,          0,          0,          0,          0,          0,              0,          0,          0,          0,          0,          0,              0,          0,          55.5248+3,    0,          0,          0,          0,              0,          0,          0,          0,          0,              0,          0,          0,          0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0],

# sm_f_in_kgm_dolo      (Flow rate of dolomite)
        [0.6304934, 0,          0,          0,              0,              0,              0,              0,          0,          0,          0,          0,          0,          0,          0,          0,          0,              0,          0.00E+00,   0,          0,          0,          0,              0,          0,          0,          0,          0,          0,          0,              0,          0,          0,          0,          0.032685746,    0.15117158, 0.145043,   0.15934302, 0.004085718,    0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0.145043,       0.10622868,     0,              0,              0,              0,              0,              0,              0,              0,              0],

# gs_f_in_ch4           (CH4 input)
        [510.04327, 1103.3473,  1017.8209,  906.6649,       764.7817,       796.7833,       818.8729,       714.08887,  195.26642,  185.77922,  183.37202,  183.65523,  182.09763,  182.52243,  183.37202,  183.23042,  922.8073,       1262.5057,  879.9025,   474.07687,  473.22726,  303.59042,  203.33763,      222.73683,  381.04565,  593.7289,   589.6225,   375.80646,  1091.7361,  1026.3169,      940.0825,   766.1977,   791.6857,   738.1609,   586.50726,      480.02405,  472.66086,  473.22726,  185.92082,      223.58643,      256.86243,      261.11044,      260.68564,      258.98642,      261.25204,      260.54404,      258.56165,      260.82724,      260.54404,      262.10162,      261.96002,      262.38483,      262.24323,      259.55283,      257.14563,      260.26083,      220.75443,      216.08163,      216.22322,      182.09763],

# ss_p_arc              (Arc power)
        [0,         49.40971,   46.85214,   55.375714,      61.421797,      61.531905,      58.308665,      58.5439,  0,          0,          0,          0,          0,          0,          0,          0,            26.266428,      63.26365,   63.89429,   66.04645,   66.00000,   66.05646,   69.15458,       67.33775, 0,          0,          0,          0,            44.82009,   52.067384,      56.28663,   54.61995,   67.85827,   70.36079,   69.83526,       64.18958,   66.34174,   72.38282,   71.73718,       65.40581,       67.793205,      75.716175,      76.9324,        76.00647,       71.93237,       67.808226,      67.30772,       69.09952,       63.67406,       63.92431,       63.84924,       63.17356,       63.93933,       64.57998,       62.67806,       61.677054,      61.96234,       61.271645,      60.575943,    0],

# sm_f_in_o2_jetbox2    (Oxygen from jetbox 1)
        [194.27522, 683.3617,   614.9689,   696.2473,       848.4673,       843.3697,       883.1593,       844.7857,   177.28322,  137.91841,  129.13922,  127.156815, 124.89121,  127.156815, 126.30721,  126.87362,  916.8601,       933.0025,   1169.1914,  1.19E+03,   1204.7329,  552.52325,  243.69363,      147.54723,  257.85364,  378.21365,  403.13525,  187.05362,  694.6897,   612.84485,      770.7289,   842.2369,   900.1513,   852.85693,  999.4129,       1166.6426,  1192.2721,  1192.2721,  342.24725,      161.14082,      141.45842,      136.50241,      127.44002,      144.57362,      138.62642,      141.60002,      128.85602,      132.82082,      151.79521,      132.96242,      136.21922,      151.08722,      132.39601,      141.31682,      123.19202,      126.87362,      143.01602,      151.22882,      143.86562,      124.04162],

# sm_f_in_o2_jetbox3    (Oxygen from jetbox 2)
        [183.79683, 694.9729,   621.3409,   740.7097,       852.0073,       848.4673,       849.6001,       846.7681,   185.21283,  164.39761,  176.43362,  163.97282,  169.77843,  169.35362,  173.74323,  161.70721,  641.87286,      876.2209,   1198.9274,  1.19E+03,   1185.0505,  743.5417,   376.37283,      807.2617,   291.97925,  389.11685,  391.38245,  179.40723,  689.4505,   624.45605,      754.5865,   853.9897,   848.6089,   852.85693,  386.85126,      1199.777,   1191.9889,  1191.9889,  392.51526,      791.2609,       1191.9889,      1931.5658,      1951.8147,      1941.9026,      1944.4514,      1937.9379,      1946.009,       1935.8138,      1952.9474,      1946.5754,      1965.833,       1940.9115,      1953.5138,      2004.7731,      1944.8762,      1940.9115,      286.59845,      168.22083,      168.36243,      171.61922],

# sm_f_in_o2_jetbox1    (Oxygen from jetbox 3)
        [207.01923, 685.76886,  621.3409,   691.0081,       851.5825,       850.5913,       860.3617,       850.7329,   166.80482,  146.55602,  149.95442,  150.09602,  148.53842,  150.37921,  148.82162,  151.22882,  740.8513,       937.10895,  1.18E+03,   1191.4226,  1191.1394,  659.4313,   1124.7289,      1942.7522,  294.38644,  388.55045,  391.24084,  191.16002,  685.2025,   622.7569,       747.08167,  852.1489,   869.28253,  849.8833,   1010.1745,      1187.5994,  1190.8562,  1189.7234,  432.02167,      1844.1986,      1965.125,       1952.8058,      1959.461,       1955.6378,      1955.213,       1945.5842,      1946.1506,      1955.7794,      1956.4874,      1957.337,       1957.337,       1939.3539,      1954.3634,      1961.5851,      1959.7443,      1952.9474,      1964.7002,      1951.6731,      1950.3987,      197.53203]])

#%%##########################################
# Matrix of inputs (Total 12)
#############################################

# Time:  0          1           2           3               4               5               6               7           8           9           10          11          12          13          14          15          16              17          18          19          20          21          22              23          24          25          26          27          28          29              30          31          32          33          34              35          36          37          38              39              40              41              42              43              44              45              46              47              48              49              50              51              52              53              54              55              56              57              58              59
# gs_f_in_kgm_oil       (volatile input)
U1 = [[0.15,      0,          0,          0,              0,              0,              0,              0,          0,          0,          0,          0,          0,          0,          0,          0,          0,              0,          0,          0,          0,          0,          0,              0,          0,          0.15,       0,          0,          0,          0,              0,          0,          0,          0,          0,              0,          0,          0,          0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0],

# sm_f_in_cao           (Flow rate of lime)
        [0.7983226, 0,          0,          0,              0,              0,              0,              0,          0.01238307, 0,          0,          0,          0,          0,          0,          0,          0,              0.27210885, 0.14912872, 0.16342874, 0.14708586, 0.08375723, 0,              0,          0,          0,          0,          0,          0,          0,              0,          0,          0,          0,          0,              0,          0,          0,          0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0],

# sm_f_in_kgm_c_charge  (Carbon charged additions)
        [1.102229,  0,          0,          0,              0,              0,              0,              0,          0,          0,          0,          0,          0,          0,          0,          0,          0,              0,          0,          0,          0,          5.00E-04,   0.029,          0.04,       0.0075,     0,          0,          0,          0,          0,              0,          0,          0,          0,          0,              0,          0,          0,          0.00472973,     0.017792793,    0.018693693,    0.02027027,     0.022972973,    0.023198199,    0.023873873,    0.012387387,    0.011711712,    0.011486487,    0.011711712,    0.013738738,    0.010135135,    0.011036036,    0.012162162,    0.012612613,    0.015990991,    0.016441442,    0.016666668,    0.017117117,    0.017117117,    0.002927928],

# sm_f_in_kgm_c_lance   (Carbon lanced additions)
        [0,         0,          0,          0,              0,              0,              0,              0,          0,          0,          0,          0,          0,          0,          0,          0,          0,              0,          0,          0,          0,          9.53E-03,   0.026535425,    0.03538057, 0,          0,          0,          0,          0,          0,              0,          0,          0,          0,          0,              0,          0,          0,          0.018370679,    0.037194956,    0.037421755,    0.04105053,     0.0460401,      0.049668875,    0.05012247,     0.023587046,    0.022453053,    0.023133447,    0.023587046,    0.028349815,    0.020638665,    0.020411866,    0.025401434,    0.025401434,    0.032658987,    0.03311258,     0.034246575,    0.03333938,     0.03515377,     0],
# gs_f_h20_electrode    (Water from sprays)
        [7.723762,  435.02954,  408.90503,  408.22354,      408.33713,      397.66016,      397.09222,      403.79373,  0,          0,          0,          0,          0,          0,          0,          0,          298.84143,      423.78464,  415.0386,   412.99408,  4.09E+02,   406.9741,   396.29715,      409.70013,  492.617,    580.3044,   578.9414,   578.9414,   393.79828,  417.76465,      418.9005,   417.4239,   411.06314,  402.99863,  394.82053,      404.47525,  415.49295,  416.7424,   417.76465,      408.90503,      395.50204,      398.56882,      408.10995,      417.65106,      419.12766,      411.85825,      397.88733,      396.29715,      409.1322,       417.87823,      420.26352,      414.4707,       402.08997,      395.1613,       402.5443,       412.53976,      419.24124,      419.6956,       415.60654,      1.9309405],

# ss_m_in_t             (2nd charge addition)
        [0,         0,          0,          0,              0,              0,              0,              0,          0,          0,          0,          0,          0,          0,          0,          0,          0,              0,          0,          0,          0,          0,          0,              0,          0,          55.5248+3,    0,          0,          0,          0,              0,          0,          0,          0,          0,              0,          0,          0,          0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0],

# sm_f_in_kgm_dolo      (Flow rate of dolomite)
        [0.6304934, 0,          0,          0,              0,              0,              0,              0,          0,          0,          0,          0,          0,          0,          0,          0,          0,              0,          0.00E+00,   0,          0,          0,          0,              0,          0,          0,          0,          0,          0,          0,              0,          0,          0,          0,          0.032685746,    0.15117158, 0.145043,   0.15934302, 0.004085718,    0,              0,              0,              0,              0,              0,              0,              0,              0,              0,              0.145043,       0.10622868,     0,              0,              0,              0,              0,              0,              0,              0,              0]]

# gs_f_in_ch4           (CH4 input)
U =     [[510.04327, 1103.3473,  1017.8209,  906.6649,       764.7817,       796.7833,       818.8729,       714.08887,  195.26642,  185.77922,  183.37202,  183.65523,  182.09763,  182.52243,  183.37202,  183.23042,  922.8073,       1262.5057,  879.9025,   474.07687,  473.22726,  303.59042,  203.33763,      222.73683,  381.04565,  593.7289,   589.6225,   375.80646,  1091.7361,  1026.3169,      940.0825,   766.1977,   791.6857,   738.1609,   586.50726,      480.02405,  472.66086,  473.22726,  185.92082,      223.58643,      256.86243,      261.11044,      260.68564,      258.98642,      261.25204,      260.54404,      258.56165,      260.82724,      260.54404,      262.10162,      261.96002,      262.38483,      262.24323,      259.55283,      257.14563,      260.26083,      220.75443,      216.08163,      216.22322,      182.09763],

# ss_p_arc              (Arc power)
        [0,         49.40971,   46.85214,   55.375714,      61.421797,      61.531905,      58.308665,      58.5439,  0,          0,          0,          0,          0,          0,          0,          0,            26.266428,      63.26365,   63.89429,   66.04645,   66.00000,   66.05646,   69.15458,       67.33775, 0,          0,          0,          0,            44.82009,   52.067384,      56.28663,   54.61995,   67.85827,   70.36079,   69.83526,       64.18958,   66.34174,   72.38282,   71.73718,       65.40581,       67.793205,      75.716175,      76.9324,        76.00647,       71.93237,       67.808226,      67.30772,       69.09952,       63.67406,       63.92431,       63.84924,       63.17356,       63.93933,       64.57998,       62.67806,       61.677054,      61.96234,       61.271645,      60.575943,    0],

# ss_p_arc              (Arc power)
       # [0,         0.04940971, 0.04685214, 0.055375714,    0.061421797,    0.061531905,    0.058308665,    0.0585439,  0,          0,          0,          0,          0,          0,          0,          0,          0.026266428,    0.06326365, 0.06389429, 0.06604645, 6.60E-02,   0.06605646, 0.06915458,     0.06733775, 0,          0,          0,          0,          0.04482009, 0.052067384,    0.05628663, 0.05461995, 0.06785827, 0.07036079, 0.06983526,     0.06418958, 0.06634174, 0.07238282, 0.07173718,     0.06540581,     0.067793205,    0.075716175,    0.0769324,      0.07600647,     0.07193237,     0.067808226,    0.06730772,     0.06909952,     0.06367406,     0.06392431,     0.06384924,     0.06317356,     0.06393933,     0.06457998,     0.06267806,     0.061677054,    0.06196234,     0.061271645,    0.060575943,    0],

# sm_f_in_o2_jetbox2    (Oxygen from jetbox 1)
        [194.27522, 683.3617,   614.9689,   696.2473,       848.4673,       843.3697,       883.1593,       844.7857,   177.28322,  137.91841,  129.13922,  127.156815, 124.89121,  127.156815, 126.30721,  126.87362,  916.8601,       933.0025,   1169.1914,  1.19E+03,   1204.7329,  552.52325,  243.69363,      147.54723,  257.85364,  378.21365,  403.13525,  187.05362,  694.6897,   612.84485,      770.7289,   842.2369,   900.1513,   852.85693,  999.4129,       1166.6426,  1192.2721,  1192.2721,  342.24725,      161.14082,      141.45842,      136.50241,      127.44002,      144.57362,      138.62642,      141.60002,      128.85602,      132.82082,      151.79521,      132.96242,      136.21922,      151.08722,      132.39601,      141.31682,      123.19202,      126.87362,      143.01602,      151.22882,      143.86562,      124.04162],

# sm_f_in_o2_jetbox3    (Oxygen from jetbox 2)
        [183.79683, 694.9729,   621.3409,   740.7097,       852.0073,       848.4673,       849.6001,       846.7681,   185.21283,  164.39761,  176.43362,  163.97282,  169.77843,  169.35362,  173.74323,  161.70721,  641.87286,      876.2209,   1198.9274,  1.19E+03,   1185.0505,  743.5417,   376.37283,      807.2617,   291.97925,  389.11685,  391.38245,  179.40723,  689.4505,   624.45605,      754.5865,   853.9897,   848.6089,   852.85693,  386.85126,      1199.777,   1191.9889,  1191.9889,  392.51526,      791.2609,       1191.9889,      1931.5658,      1951.8147,      1941.9026,      1944.4514,      1937.9379,      1946.009,       1935.8138,      1952.9474,      1946.5754,      1965.833,       1940.9115,      1953.5138,      2004.7731,      1944.8762,      1940.9115,      286.59845,      168.22083,      168.36243,      171.61922],

# sm_f_in_o2_jetbox1    (Oxygen from jetbox 3)
        [207.01923, 685.76886,  621.3409,   691.0081,       851.5825,       850.5913,       860.3617,       850.7329,   166.80482,  146.55602,  149.95442,  150.09602,  148.53842,  150.37921,  148.82162,  151.22882,  740.8513,       937.10895,  1.18E+03,   1191.4226,  1191.1394,  659.4313,   1124.7289,      1942.7522,  294.38644,  388.55045,  391.24084,  191.16002,  685.2025,   622.7569,       747.08167,  852.1489,   869.28253,  849.8833,   1010.1745,      1187.5994,  1190.8562,  1189.7234,  432.02167,      1844.1986,      1965.125,       1952.8058,      1959.461,       1955.6378,      1955.213,       1945.5842,      1946.1506,      1955.7794,      1956.4874,      1957.337,       1957.337,       1939.3539,      1954.3634,      1961.5851,      1959.7443,      1952.9474,      1964.7002,      1951.6731,      1950.3987,      197.53203]]

# We transpose the U matrix for later use
U_dm2 = [[row[i11] for row in U1] for i11 in range(60)]
U_dm1 = [[row[i111] for row in U] for i111 in range(60)]

# Matrix of the model noises (type of input)
w_dm = [[0]*len(dae.x)]*60

# Eliminate the algebraic variables (not all will be eliminated !)
dae.eliminate_alg()

dae_x = []
for i in range(len(dae.x)):
    dae_x = vertcat(dae_x,dae.x[i])

dae_z = []
for i_ in range(len(dae.z)):
    dae_z = vertcat(dae_z,dae.z[i_])
    
dae_p = []
for i1_ in range(len(dae.p)):
    dae_p = vertcat(dae_p,dae.p[i1_])    
for i2_ in range(len(dae.u)):
    dae_p = vertcat(dae_p,dae.u[i2_]) 
    
dae_ode = []
for i in range(len(dae.ode)):
    dae_ode = vertcat(dae_ode,dae.ode[i])    

dae_alg = []
for i_ in range(len(dae.alg)):
    dae_alg = vertcat(dae_alg,dae.alg[i_])    

#f = dae.create('f',['x','z','u','p'],['ode','alg'])
dae1 = {'x':dae_x, 'z':dae_z, 'p':dae_p, 'ode':dae_ode, 'alg':dae_alg}#, 'quad':quad}

Z_reduced = []
for i3_ in range(len(dae.alg)):
    Z_reduced = Z_reduced + dae.guess(dae.z[i3_]) 

#%%###################################################
#          Theretical optima solve starts here
######################################################



#%%###################################################
#                  MHE starts here
######################################################

## Linearization calculations
# d = Ex + Fz
# 0 = Gx + Hz
# x_dot = Ix + Jz
# z = -inv(H)*G*x
# xdot = (I - J*inv(H)*G)x
# d = (E - F*inv(H)*G)x

# E,F,G,H,I,J
f_c = dae.create('f_c', ['x','z','u','p'], ['jac_ddef_x','jac_ddef_z',
                 'jac_alg_x','jac_alg_z','jac_ode_x','jac_ode_z'])

## Define a measurement function
f_y = dae.create('f_y',['x','z','u','p'],['ddef'])

# Total MHE solves (including batch MHE solves)
n_mhe = 60

# MHE horizon
mhe_horizon = 6 #(total 7 measurements in the horizon)

# start of MHE
n_ini = 0

# Get the actual initial states
x_ini_start_mhe = X0
z_ini_start_mhe = Z_reduced

# Get the number of states and algebraic variables
Nstates = len(x_ini_start_mhe)
Nalgvars = len(z_ini_start_mhe)
Ninputs = 12

# Column vectors to be used later in C
if disturb_add == 1:
    my_c_t_roof = horzcat(DM.zeros(1, Nstates - 3), 1, 0, 0)
    my_c_t_wall = horzcat(DM.zeros(1, Nstates - 3), 0, 1, 0)
    my_c_mm_t = horzcat(DM.zeros(1, Nstates - 10), 1,
                        DM.zeros(1, Nstates - 20))
else:                        
    my_c_t_roof = horzcat(DM.zeros(1, Nstates - 2), 1, 0)
    my_c_t_wall = horzcat(DM.zeros(1, Nstates - 2), 0, 1)
    my_c_mm_t = horzcat(DM.zeros(1, Nstates - 9), 1,
                        DM.zeros(1, Nstates - 20))                   

# Measurement structure
# Time:        0    1   2   3   4   5   6   7   8   9   10  11  12  13  14  15    
meas_struct = [6,   6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,
# Time:        16   17  18  19  20  21  22  23  24  25  26  27  28  29  30  31               
               6,   6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,
# Time:        32   33  34  35  36  37  38  39  40  41  42  43  44  45  46  47       
               #6,   6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  13,  6,  6,  6,  8,
               6,   6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  13,  6,  6,  6,  8,
# Time:        48   49  50  51  52  53  54  55  56  57  58  59  60       
               6,   6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6]
                                   
#%% Step 1: Add noise on the initial states  
sigma_x0 = DM.zeros(Nstates,1)
#sigma_x1 = DM.zeros(Nstates,1)
for i in range(Nstates):
    # make sure you do not add more than 5 and less than 0.01
    #sigma_x0[i] = max(0.0001,min(abs(0.01*x_ini_start_mhe[i]),5))
    ## sigma_x0[i] = max(0.0001,min(abs(0.01*x_ini_start_mhe[i]),10))  # Final
    sigma_x0[i] = max(0.0001,min(abs(0.01*x_ini_start_mhe[i]),10))  # Final
    #sigma_x1[i] = max(0.00001,min(abs(0.10*x_ini_start_mhe[i]),5))  # Final
    
P = DM.eye(Nstates)
for i in range(Nstates):
    P[i,i] = sigma_x0[i]**2

S = linalg.inv(P) # Arrival cost is the inverse of the initial covariance    

# Get the noise added initial states 
NP.random.seed(0) 
x_ini_start_mhe_wn = x_ini_start_mhe + sigma_x0*NP.random.randn(Nstates,1)

# Reset known stats values
x_ini_start_mhe_wn[8:11] = 0
x_ini_start_mhe_wn[24] = 0
if disturb_add ==1:
    x_ini_start_mhe_wn[Nstates-1] = 0
 
# Build the state matrices
x_mat_act = []
x_mat_act = x_ini_start_mhe

z_mat_act = []
z_mat_act = z_ini_start_mhe

x_mat_est = []
x_mat_est = x_ini_start_mhe_wn

#%%####################################################################### 
# Solve the theoretical and nominal optimization problems
#######################################################################

# Price change time 
n_change = 40

# Profit function coefficients
# Time invariant cost for CH4
u7_coef = 1.5*0.213650/60 # CH4 inupt

#u8_coef = 0.0585*1000/60  # Arc Power
#u8_coef = 0.0350*1000/60  # Arc Power
#u8_coef = []  # 35 $/MWh = 35/60 $/MWmin = 35*0.001/60 

# We are having a change in price at time t = 30 min, units: $/MWh
u8_coef1_half = [308.24]*n_change
# u8_coef1 = u8_coef1_half[30] + [15]*(n_mhe-n_change), 200 afterwards # 75-76
# u8_coef1 = u8_coef1_half[308.24] + [10.96]*(n_mhe-n_change), 190.48 afterwards
# 30 first, 41.85 predicted, actual 20.53
u8_coef1 = u8_coef1_half + [10.96]*(n_mhe-n_change)
u8_coef = [x1_ * (1/60) for x1_ in u8_coef1]
# Updated price at t=30 min
#u8_coef1_updated = [38.5]*(n_mhe-n_change)
u8_coef1_updated = [190.48]*(n_mhe-n_change)
# Real price profile
u8_coef1_updated1 = [190.48]*(n_mhe-n_change)
u8_coef21 = u8_coef1_half + u8_coef1_updated1
u8_coef2 = [x1_ * (1/60) for x1_ in u8_coef21]

# Time invariant cost for O2
u9_coef = 1.5*0.07875/60  # For u10 and u11 as well
#u9_coef = 1.5*0.07875/60  # For u10 and u11 as well

xk15o_coef = 0.045*10# 0.359472*1 can be increased if required 
xk15i_coef = 55.85
xk15g_coef = 1 #0.1

# Replace the last element of initial x with 0
o_x_ini_start_mhe_wn = x_ini_start_mhe_wn

# EMPC solvetime container
solve_time_opt = []
    
for o_solve in range(2):  # 0: Theoretical, 1: Nominal
    # Define mode:- 0: feasibility/ 1:optimization
    o_mode  = 1
    
    # Specify the time where to start optimization
    o_n_ini = 0  # 2 working fine
    
    # Accordingly initialize the states and algebraic variables		
    # if o_solve == 0:
    o_X_ini = DM(x_ini_start_mhe)     # X0
        #o_X_ini_be = DM(x_ini_start_mhe)  # X0
    # else:
    o_X_ini1 = o_x_ini_start_mhe_wn   
        #o_X_ini_be = o_x_ini_start_mhe_wn 
    
    # We are trying to replicate Case Study 1 of Richard's Thesis 
    # 1) Same initial values
    # 2) Mismatch in k_p or k_dm					
    
    o_Z_ini = Z_reduced                                          
    
    # No of control intervals
    o_N = 60 - o_n_ini

    #%%##########################################
    # Setup the optimization here
    #############################################
    
    # No of steps per control interval
    o_steps = 4 #(4 works fine),(1,2,3,5,6locally infeasible),(7:max iterations):'limited-memory'ipopt:N=60 
    # 21/4 :( 25/25/5 :( (after that only 4 & 8)
            
    # step size
    o_h = 1/o_steps
    
    # output grid
    time_grid = linspace(0, 1, o_steps+1)
    
    # setup the idas integrator
    I_full = integrator('integrator', 'idas', dae1, {'tf':1, 'calc_ic':True,
                                               'grid': time_grid, 'output_t0':True})
    
    # Setup the backward euler integration
    f = dae.create('f',['x','z','p','u'],['ode','alg'])
    
    # Start with an empty NLP
    w=[]
    w0 = []
    o_lbw = []
    o_ubw = []
    J = 0
    g=[]
    lbg = []
    ubg = []
    
    # "Lift" initial conditions
    o_X_ini_def = SX.sym('o_X_ini_def', Nstates)
    w += [o_X_ini_def]
    if o_solve == 0:
	    o_lbw += o_X_ini.nonzeros()
	    o_ubw += o_X_ini.nonzeros()
	    w0 += o_X_ini.nonzeros()				
    else:
	    o_lbw += o_X_ini1.nonzeros()
	    o_ubw += o_X_ini1.nonzeros()
	    w0 += o_X_ini1.nonzeros()
    
    o_Xk = o_X_ini_def
    
    # Define contrainer for exact state values obtained from IDAS
    if o_solve == 0: 
	    o_X_exact = o_X_ini
	    o_Z_exact = o_Z_ini
    else: 
	    o_X_exact = o_X_ini1
	    o_Z_exact = o_Z_ini
					
    #%%############################################################
    # Back euler loop starts here
    ############################################################### 
    for o_j1 in range(o_N):
        
        # Initialization of backeuler simulation
        Ik1 = I_full(x0=o_X_ini, 
                      p=vertcat(par_val_mismatch, 
                                U_dm2[o_j1+o_n_ini] + U_dm1[o_j1+o_n_ini] + w_dm[o_j1+o_n_ini]), 
                      z0 = o_Z_ini)
																					
        o_X_ini = Ik1['xf'][:,-1]
        o_Z_ini = Ik1['zf'][:,-1]
        o_X_exact = horzcat(o_X_exact, o_X_ini)
        o_Z_exact = horzcat(o_Z_exact, o_Z_ini)
        
        # NLP part
        # New NLP variable for the control
        o_Uk = SX.sym('o_U_' + str(o_j1), Nstates + Ninputs)
        w   += [o_Uk]
        
        # Fix inputs according to the feasibility/optimization
        if o_mode == 0:
            # Inputs fixed as we want to solve for feasibility only
            o_lbw += U_dm2[o_j1+o_n_ini] + U_dm1[o_j1+o_n_ini]
            o_ubw += U_dm2[o_j1+o_n_ini] + U_dm1[o_j1+o_n_ini]
        else:    
            o_lbw += [x1_ * 1 for x1_ in U_dm2[o_j1+o_n_ini]]
            o_lbw += [0.9*U_dm1[o_j1+o_n_ini][0]]
            o_lbw += [0.9*U_dm1[o_j1+o_n_ini][1]]#ss_p_arc 0.8
            o_lbw += [x2_ * 0.9 for x2_ in U_dm1[o_j1+o_n_ini][2:]]
        
            o_ubw += [x3_ * 1 for x3_ in U_dm2[o_j1+o_n_ini]]
            o_ubw += [1.1*U_dm1[o_j1+o_n_ini][0]]  
            o_ubw += [1.1*U_dm1[o_j1+o_n_ini][1]] #ss_p_arc 1.2
            o_ubw += [x4_ * 1.1 for x4_ in U_dm1[o_j1+o_n_ini][2:]]
            
            # Add contribution of inputs to cost function
            J = J  + u7_coef*o_Uk[7] + (u8_coef[o_j1]*o_Uk[8] + u9_coef*(o_Uk[9]+o_Uk[10]+o_Uk[11]))
            #J = J  + 1.5*0.213650*Uk[7]/60 + (0.0585*1000*Uk[8]/60 + 1.5*0.07875*(Uk[9]+Uk[10]+Uk[11])/60)

        # Provide bounds for model noises 
        o_lbw += [0]*Nstates
        o_ubw += [0]*Nstates
            
        w0  += U_dm2[o_j1+o_n_ini] + U_dm1[o_j1+o_n_ini] + w_dm[o_j1+o_n_ini]
        
        for o_j in range(o_steps):
    
            o_X_ini_be = Ik1['xf'][:,o_j+1]
            o_Z_ini_be = Ik1['zf'][:,o_j+1]
            
            # NLP part
            # State at the euler points
            o_Xk_end = SX.sym('o_X_' + str(o_j1+1) + '_' + str(o_j), Nstates)
            o_Zk_end = SX.sym('o_Zz_' + str(o_j1+1) + '_' + str(o_j), Nalgvars)
            
            w   += [o_Xk_end]
    
            for i5 in range(Nstates):
                if i5 in range(Nstates-1):#[26,27,7]:#,7,8,10]:
                    if abs(o_X_ini_be[i5]) > 0 and abs(o_X_ini_be[i5]) < 5:
                        if o_X_ini_be[i5] < 0:
                            o_lbw += [-10]
                        else:
                            o_lbw += [0]
                    else:
                        if o_X_ini_be[i5] <= 0:
                            o_lbw += (1.5*o_X_ini_be[i5]).nonzeros()
                        else:    
                            o_lbw += (0.5*o_X_ini_be[i5]).nonzeros()
                else:  # 0 for disturbance state   
                    # lbw+= [X_ini_be.nonzeros()[i5]] 
                    o_lbw += [0]
               
            for i4 in range(Nstates):
                if i4 in range(Nstates-1):#[26,27,7]:#,7,8,10]:            
                    if abs(o_X_ini_be[i4]) > 0 and abs(o_X_ini_be[i4]) < 5:
                        if o_X_ini_be[i4] < 0:
                            o_ubw += [0]
                        else:
                            o_ubw += [10]
                    else:
                        if o_X_ini_be[i4] < 0:
                            o_ubw += (0.5*o_X_ini_be[i4]).nonzeros()
                        else:    
                            o_ubw += (1.5*o_X_ini_be[i4]+0.0001).nonzeros()  
                else:  # 0 for disturbance state
                    # ubw += [X_ini_be.nonzeros()[i4]] 
                    o_ubw += [0]
                    
            w0 += [k for k in o_X_ini_be.nonzeros()]
            
            w   += [o_Zk_end]
        
            for i3 in range(Nalgvars):
                if i3 in range(14)+[15,17,18,20,22,24,26]+range(28,43)+range(52,108)+[110,113,115,119,120]: 
                    if abs(o_Z_ini_be[i3]) > 0 and abs(o_Z_ini_be[i3]) < 5:
                        if o_Z_ini_be[i3] < 0:
                            o_lbw += [-10]
                        else:
                            o_lbw += [0]
                    else:
                        if o_Z_ini_be[i3] <= 0:
                            o_lbw += (1.5*o_Z_ini_be[i3]).nonzeros()
                        else:    
                            o_lbw += (0.5*o_Z_ini_be[i3]).nonzeros()
                #else:      #  0     1   2   3    4    5   6   7    8   9   10  11   12  13   14   15     
                #    if i3 in [14*, 51, 16, 44*, 43*, 19, 45, 21*, 46, 23*, 47, 25*, 48, 27*, 49*, 50*]
                else:      #  0     1   2   3    4    5   6   7    8   9   10  11   12  13   14   15     
                    if i3 in [51, 16, 19, 45, 46, 47, 48]:
                        o_lbw += (1.2*o_Z_ini_be[i3]).nonzeros()
                    else:
                        o_lbw += (1.5*o_Z_ini_be[i3]).nonzeros()
                        #lbw += [-1e10]
    
            for i2 in range(Nalgvars):
                if i2 in range(14)+[15,17,18,20,22,24,26]+range(28,43)+range(52,108)+[110,113,115,119,120]: 
                    if abs(o_Z_ini_be[i2]) > 0 and abs(o_Z_ini_be[i2]) < 5:
                        if o_Z_ini_be[i2] < 0:
                            o_ubw += [0]
                        else:
                            o_ubw += [10]
                    else:
                        if o_Z_ini_be[i2] < 0:
                            o_ubw += (0.5*o_Z_ini_be[i2]).nonzeros()
                        else:    
                            o_ubw += (1.5*o_Z_ini_be[i2]+0.0001).nonzeros()
                else:      #  0     1   2   3    4    5   6   7    8   9   10  11   12  13   14   15     
                    if i2 in [51, 16, 19, 45, 46, 47, 48]:
                        o_ubw += (0.8*o_Z_ini_be[i2]).nonzeros()
                    else:
                        #ubw += (0.9*Z_ini_be[i2]).nonzeros()                    
                        #ubw += [0.001] # [0.0001] upto 47      
                        o_ubw += [0.0001]# upto 47                       
                        
            w0  += [k for k in o_Z_ini_be.nonzeros()]
        
            # Euler equations
            #f = Function('f', [x_all, z_all, params], [ode, alg])
            if o_solve == 0:  #  Theoretical
                o_fj, o_zj = f(o_Xk_end, o_Zk_end, par_val, o_Uk) 
            else:             #  Nominal
                o_fj, o_zj = f(o_Xk_end, o_Zk_end, par_val_mismatch, o_Uk)
                
            g += [o_Xk[0:(Nstates-1)] - o_Xk_end[0:(Nstates-1)] + o_fj[0:(Nstates-1)]*o_h]
            
            lbg += [0]*(Nstates-1)
            ubg += [0]*(Nstates-1)
            
            g += [o_zj]
            lbg += [0]*Nalgvars
            ubg += [0]*Nalgvars        
        
            # Add equality constraint
            o_Xk = o_Xk_end
    
    if o_mode == 1: # if mode is optimization then add objective contribution
        #J = J - 0.359472*Xk_end[12]

        J = J - xk15o_coef*((o_Xk_end[15]*xk15i_coef*xk15g_coef) + (((o_Zk_end[42])*exp(o_Zk_end[44]))*xk15i_coef))

        #J = J - 0.359472*((o_Xk_end[15]*55.85*0.1) + (((o_Zk_end[42])*exp(o_Zk_end[44]))*55.85))
    #%%#################################################
    # Solve the NLP
    ###################################################    
    
    # NLP solver options
    opts_nlp = {}
    opts_nlp["ipopt.linear_solver"] = 'ma27'
    #opts_nlp["ipopt.linear_system_scaling"] = 'slack-based'
    #opts_nlp["ipopt.alpha_for_y"]='primal-and-full'
    #opts_nlp["ipopt.obj_scaling_factor"] = 1e-3
    opts_nlp["ipopt.hessian_approximation"] = 'limited-memory'
    if o_solve == 1:  # warm start the solve
        opts_nlp["ipopt.warm_start_init_point"] = 'yes'
        opts_nlp["ipopt.warm_start_bound_push"] = 1e-6
        opts_nlp["ipopt.warm_start_mult_bound_push"] = 1e-6
        opts_nlp["ipopt.mu_init"] = 1e-6  
        
    prob = {'f': J, 'x': vertcat(*w), 'g': vertcat(*g)}
    solver = nlpsol('solver', 'ipopt', prob, opts_nlp);
    #solver = nlpsol("solver", "sqpmethod", prob, {"qpsol":"qpoases",#})
    #                                            "hessian_approximation":"limited-memory",
    #                                             "print_header":True,"max_iter":100,
    #                                             "verbose":True})#, opts)                                            "hessian_approximation":"limited-memory"})
    
    
    # Start clock to measure simulation time
    startTime = time.time()
    
    # Solve the NLP
    if o_solve == 1:  # warm start the solve
        opt_sol_nom = solver(x0=opt_sol_theo['x'], lbx=o_lbw, ubx=o_ubw, lbg=lbg, 
        #             ubg=ubg, lam_x0=opt_sol_theo['lam_x'],
        #             lam_g0 = opt_sol_theo['lam_g'])
        #opt_sol_nom = solver(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, 
                     ubg=ubg, lam_x0=opt_sol_theo['lam_x'],
                     lam_g0 = opt_sol_theo['lam_g'])
    else:
        opt_sol_theo = solver(x0=w0, lbx=o_lbw, ubx=o_ubw, lbg=lbg, ubg=ubg)                 
    
    # Stop the clock and read in the total time taken                    
    elapsedTime = time.time() - startTime
    
    if o_solve == 1:
        # Build the solve time array
        solve_time_opt = append(solve_time_opt, elapsedTime)			
        o_w_opt_nom = opt_sol_nom['x'].full().flatten()
    else:   
        o_w_opt_theo = opt_sol_theo['x'].full().flatten()

# Extract the optimal inputs from both the solutions
# Repeating part
o_w_opt_nom1 = o_w_opt_nom[Nstates:]
o_w_opt_theo1 = o_w_opt_theo[Nstates:]

opt_u_nom = []
opt_u_theo = []
                     
for iw in range(o_N):
    # Select the first set of values for a finite element
    o_w_opt_nom2 = DM(o_w_opt_nom1[(Ninputs+(Nstates)*(o_steps+1)+(Nalgvars)*o_steps)*iw:\
                       (Ninputs+(Nstates)*(o_steps+1)+(Nalgvars)*o_steps)*(iw+1)])
            
    # Now extract the input values (model noise included :| )               
    opt_u_nom = horzcat(opt_u_nom, o_w_opt_nom2[0:Ninputs+Nstates])

    # Select the first set of values for a finite element
    o_w_opt_theo2 = DM(o_w_opt_theo1[(Ninputs+(Nstates)*(o_steps+1)+(Nalgvars)*o_steps)*iw:\
                       (Ninputs+(Nstates)*(o_steps+1)+(Nalgvars)*o_steps)*(iw+1)])
            
    # Now extract the input values (model noise included :| )               
    opt_u_theo = horzcat(opt_u_theo, o_w_opt_theo2[0:Ninputs+Nstates])            

# Save the nominal inputs for plotting
opt_u_nom_ol = opt_u_nom

# Implement the nominal solution on the real plant and calculate the profit
o_xstart = x_ini_start_mhe
o_zstart = Z_reduced
J_nom_op = 0

for opt_i in range(o_N):
    np_IK1 = I_full(x0=o_xstart, p=vertcat(par_val, opt_u_nom[:,opt_i]), 
                    z0 = o_zstart)
    o_xstart = np_IK1['xf'][:,-1] 
    o_zstart = np_IK1['zf'][:,-1]                      
                          
    J_nom_op = J_nom_op  + u7_coef*opt_u_nom[7,opt_i] + \
               (u8_coef[opt_i]*opt_u_nom[8,opt_i] + u9_coef*(opt_u_nom[9,opt_i]+opt_u_nom[10,opt_i]+opt_u_nom[11,opt_i]))
    
    if opt_i == o_N-1: 
        J_nom_op = J_nom_op - xk15o_coef*((np_IK1['xf'][15,-1]*xk15i_coef*xk15g_coef) + (((np_IK1['zf'][42,-1])*exp(np_IK1['zf'][44,-1]))*xk15i_coef))
                
#%% Covariance matrix for model noise
if disturb_add == 1:
    Q_array = [0.08,  # sm_b_0 (0-7)
               1,#0.7,#0.2,  # sm_b_1 (0-110)
               0.6,#0.2,#0.1,  # sm_b_2 (0-100)
               0.1,#0.01,#recent graph,#0.008,  # sm_b_3 (0.05-0.08)
               0.1,  # sm_b_4 (0-12)
               0.1,  # sm_b_5 (0-15)
               0.1,  # sm_b_6 (0-10)
               0.4,#0.3,#0.1,  # sm_b_7 (0-100)
               0.03,  # sm_m_cao (0-0.6)
               0.03,  # sm_m_c_float (0-0.8)
               0.06,  # sm_m_dol (0-0.45)
               0.4,#0.2,#0.1,  # sm_e (-70000-0)
               0.01,#0.1,  # sm_m_ss (0-55)
               0.1,#0.1,#working base 0.1,  # ss_t (0-1900)
               0.1,  # mm_n_0 (0-6)
               3,#1,#0.1,  # mm_n_1 (30-1900)
               0.1,  # mm_n_2 (0.1-1.5)
               0.7,#0.4,#0.2,#0.1,  # mm_n_3 (0-90)
               0.1,  # mm_n_4 (0-70)
               3,#1,#0.1,  # mm_t (0-2000)
               0.08,  # gs_b_0 (0-12) 
               0.1,  # gs_b_1 (0-35)
               0.1,  # gs_b_2 (0-35)
               0.1,  # gs_b_3 (0-40)
               0.01,  # gs_n_oil_gas (0-0.7)
               4,#2,#0.1,  # gs_t (300-2000)
               0.1,  # rd_t_roof (298-314)
               0.1,  # rd_t_wall (298-308)
               0.2]#0.05]#0.01]  # d_sm_b_3 
else:
    Q_array = [0.08,  # sm_b_0 (0-7)
               1,#0.7,#0.2,  # sm_b_1 (0-110)
               0.6,#0.2,#0.1,  # sm_b_2 (0-100)
               0.1,#0.01,#recent graph,#0.008,  # sm_b_3 (0.05-0.08)
               0.1,  # sm_b_4 (0-12)
               0.1,  # sm_b_5 (0-15)
               0.1,  # sm_b_6 (0-10)
               0.4,#0.3,#0.1,  # sm_b_7 (0-100)
               0.03,  # sm_m_cao (0-0.6)
               0.03,  # sm_m_c_float (0-0.8)
               0.06,  # sm_m_dol (0-0.45)
               0.4,#0.2,#0.1,  # sm_e (-70000-0)
               0.1,  # sm_m_ss (0-55)
               0.1,#0.1,#working base 0.1,  # ss_t (0-1900)
               0.1,  # mm_n_0 (0-6)
               3,#1,#0.1,  # mm_n_1 (30-1900)
               0.1,  # mm_n_2 (0.1-1.5)
               0.7,#0.4,#0.2,#0.1,  # mm_n_3 (0-90)
               0.1,  # mm_n_4 (0-70)
               3,#1,#0.1,  # mm_t (0-2000)
               0.08,  # gs_b_0 (0-12) 
               0.1,  # gs_b_1 (0-35)
               0.1,  # gs_b_2 (0-35)
               0.1,  # gs_b_3 (0-40)
               0.01,  # gs_n_oil_gas (0-0.7)
               4,#2,#0.1,  # gs_t (300-2000)
               0.1,  # rd_t_roof (298-314)
               0.1]  # rd_t_wall (298-308)
               #0.05]#0.01]  # d_sm_b_3 
           
sigma_w = DM.zeros(Nstates,1)
for i in range(Nstates):
    sigma_w[i] = Q_array[i]#0.1

Q1 = DM.eye(Nstates)
for i in range(Nstates):
    Q1[i,i] = sigma_w[i]**2
    
Q = linalg.inv(Q1)

# Bounds for the model noises
if disturb_add == 1:
    w_ub = DM([0.5,   # sm_b_0 (0-7)
               1,     # sm_b_1 (0-110)
               1,     # sm_b_2 (0-100)
               0.5,#0.01,#0.008,  # sm_b_3 (0.05-0.08) 
               0.5,   # sm_b_4 (0-12) 
               0.5,   # sm_b_5 (0-15) 
               0.4,   # sm_b_6 (0-10) 
               1,     # sm_b_7 (0-100)
               0.08,#0.3,  # sm_m_cao (0-0.6)
               0.3,  # sm_m_c_float (0-0.8)
               0.6,  # sm_m_dol (0-0.45)
               1,     # sm_e (-70000-0)
               1,     # sm_m_ss (0-55)
               1,     # ss_t (0-1900)
               0.3,   # mm_n_0 (0-6)
               1,     # mm_n_1 (30-1900)
               0.8,  # mm_n_2 (0.1-1.5)
               1,     # mm_n_3 (0-90)
               1,     # mm_n_4 (0-70)
               1,     # mm_t (0-2000)
               0.5,   # gs_b_0 (0-12) 
               0.7,   # gs_b_1 (0-35)
               0.7,   # gs_b_2 (0-35)
               0.7,   # gs_b_3 (0-40)
               0.5,  # gs_n_oil_gas (0-0.7)
               1,     # gs_t (300-2000)
               1,     # rd_t_roof (298-314)
               1,    # rd_t_wall (298-308)
               0.5])   # d_sm_b_3
else:  
    w_ub = DM([0.5,   # sm_b_0 (0-7)
               1,     # sm_b_1 (0-110)
               1,     # sm_b_2 (0-100)
               #0.01,#0.008,  # sm_b_3 (0.05-0.08) 
               0.5,#0.008,  # sm_b_3 (0.05-0.08)                
               0.5,   # sm_b_4 (0-12) 
               0.5,   # sm_b_5 (0-15) 
               0.4,   # sm_b_6 (0-10) 
               1,     # sm_b_7 (0-100)
               0.08,#0.3,  # sm_m_cao (0-0.6)
               0.3,  # sm_m_c_float (0-0.8)
               0.6,  # sm_m_dol (0-0.45)
               1,     # sm_e (-70000-0)
               1,     # sm_m_ss (0-55)
               1,     # ss_t (0-1900)
               0.3,   # mm_n_0 (0-6)
               1,     # mm_n_1 (30-1900)
               0.8,  # mm_n_2 (0.1-1.5)
               1,     # mm_n_3 (0-90)
               1,     # mm_n_4 (0-70)
               1,     # mm_t (0-2000)
               0.5,   # gs_b_0 (0-12) 
               0.7,   # gs_b_1 (0-35)
               0.7,   # gs_b_2 (0-35)
               0.7,   # gs_b_3 (0-40)
               0.5,  # gs_n_oil_gas (0-0.7)
               1,     # gs_t (300-2000)
               1,     # rd_t_roof (298-314)
               1])#,    # rd_t_wall (298-308)
           
w_lb = -w_ub           
    
#%% Step 2: Add noise on the initial measurements

# Build the measurement weighing matrices
"""
sigma_y_full = DM([0.01, 0.01, 0.01, 0.01, 3, 3,
                   5, 0.01,
                   0.01, 0.01, 0.01, 0.01, 0.01]) # all measurments
"""                   

sigma_y_full = DM([0.01, 0.01, 0.01, 0.01, 3, 3,
                   5, 0.01,#0.01
                   #0.5, 0.5, 0.5, 0.5, 0.5]) # all measurments                   
                   0.1, 0.1, 0.1, 0.1, 0.1]) # all measurments
                   #0.01, 0.01, 0.01, 0.01, 0.01]) # all measurments
                   
sigma_y_half = sigma_y_full[0:8] # continuous  + mm measurments
sigma_y = sigma_y_full[0:6] # continuous only

# all measurments
R1_full = DM.eye(sigma_y_full.shape[0])
for i in range(sigma_y_full.shape[0]):
    R1_full[i,i] = sigma_y_full[i]**2
 
R_full = linalg.inv(R1_full)

# continuous  + mm measurments
R1_half = DM.eye(sigma_y_half.shape[0])
for i in range(sigma_y_half.shape[0]):
    R1_half[i,i] = sigma_y_half[i]**2
 
R_half = linalg.inv(R1_half)   

# continuous only
R1 = DM.eye(sigma_y.shape[0])
for i in range(sigma_y.shape[0]):
    R1[i,i] = sigma_y[i]**2
 
R = linalg.inv(R1)

# Define the integrators
I_first = integrator('integrator', 'idas', dae1, {'tf':1, 'calc_ic':True,
                                           #'linear_solver_type':'user_defined',
                                           'linear_solver':'csparse',
                                           'output_t0':True})

I_next = integrator('integrator', 'idas', dae1, {'tf':1, 'calc_ic':True,
                                           #'linear_solver_type':'user_defined',
							              'abstol':1e-6, 'reltol':1e-6, # earlier was 8 																																											
                                           'linear_solver':'csparse',
                                           'output_t0':False})                                           

# Simulate the model one time step and get the first measurment 
Ik_first = I_first(x0 = x_ini_start_mhe, 
                   p = vertcat(par_val, opt_u_nom[:,0]),
                   z0 = z_ini_start_mhe)

# Get all the dependent variable values
d_all = f_y(Ik_first['xf'][:,0], Ik_first['zf'][:,0], opt_u_nom[:,0], 
            par_val)

# Get the measurement count at time t=0
Nmeas = meas_struct[0] 
                  
# Choose to get the measurement vector and add noise accordingly
y_first = vertcat(d_all[318:321], d_all[322], Ik_first['xf'][26:28,0],# continuous measurements
                  Ik_first['xf'][19,0], d_all[396], # molten metal measurements 
                  d_all[196], d_all[201], d_all[203], d_all[205:207]) #slag measurements          

# Add noise
NP.random.seed(1)
y_first_wn = y_first + sigma_y_full*NP.random.randn(y_first.shape[0],1)                                    

# Build the measurement matrices and add the first one
y_mat = []
y_mat = y_first

y_mat_wn = []
y_mat_wn = y_first_wn

# Select measurments according to measurement structure at initial time   
if Nmeas == 8:                  
    y_first_wn = y_first_wn[0:8] 
    
if Nmeas == 6:                  
    y_first_wn = y_first_wn[0:6]

# Build the solve time storage
solve_time_mhe = []
solve_time_pmhe = []
solve_time_fmhe = []   

solve_time_pmpc = []
solve_time_fmpc = [] 

# Build the ipopt iteration count storage
solve_iter_mhe = []        

#%%#####################################################
#           Setup the MHE optimization here
########################################################

# No of steps per control interval
steps = 7 #5(4 works fine),(1,2,3,5,6locally infeasible),(7:max iterations):'limited-memory'ipopt:N=60 
# 21/4 :( 25/25/5 :( (after that only 4 & 8)

# Step size
h = 1/steps

# Output grid
time_grid = linspace(0, 1, steps+1)

# setup the idas integrator
I_full = integrator('integrator', 'idas', dae1, {'tf':1, 'calc_ic':True,
                                           #'linear_solver_type':'user_defined',
                                           'abstol':1e-8, 'reltol':1e-8,              
                                           'linear_solver':'csparse',
                                           'grid': time_grid, 'output_t0':True})

# Setup the backward euler integration
f = dae.create('f',['x','z','p','u'],['ode','alg'])

# Optimization initialization 
X_ini = DM(x_ini_start_mhe_wn)#_wn
Z_ini = DM(z_ini_start_mhe) # no noise added
X_ini_be = x_ini_start_mhe_wn
#X_ini_be = x_ini_start_mhe #no noise added
Z_ini_be = DM(z_ini_start_mhe) # no noise added

# Replace the last element of initial x with 0
# X_ini[8:11] = 0
# X_ini[24] = 0
# if disturb_add ==1:
#     X_ini[Nstates-1] = 0

# X_ini_be[8:11] = 0
# X_ini_be[24] = 0
# if disturb_add ==1:
#     X_ini_be[Nstates-1] = 0
    
###################
#"""
## Finding consistent guesses
f11 = Function('f11', [dae_z, dae_x, dae_p], [dae_alg, dae_z, dae_ode])

opts_test1 = {}
opts_test1["abstol"] = 1e-14
opts_test1["linear_solver"] = "csparse"
opts_test1["print_iteration"] = False

s11 = rootfinder("s1", "newton", f11, opts_test1)
x_1 = s11(Z_ini_be, X_ini_be, 
          vertcat(par_val_mismatch, 
                  opt_u_nom[:,n_ini]))

# Assign the new consistent alg var values 
Z_ini = x_1[0]
Z_ini_be = x_1[0]
#"""
###################

# Optimization related declarations
# opt_call = [10, 20, 30, 40, 50]
opt_call = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
			19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34,
			35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
			50, 51, 52, 53, 54, 55, 56, 57, 58, 59]
# diff = [opt_call[0], opt_call[1]-opt_call[0], opt_call[2]-opt_call[1],
#        opt_call[3]-opt_call[2], opt_call[4]-opt_call[3]]
diff = 1
count = 0  
J_cl = 0
J_cl_arc_power = 0 
J_cl_other = 0   

#%%###########################################################################
#                                MHE loop 
##############################################################################
for i in range(n_mhe-1): # 27 was fine, range(11) range(n_mhe): 26,35,36,37,38,39,40,41,42,43,44: max iteration but fast converging

    #%% Measurement build up 
    # Simulate the plant one step 
    Ik_next = I_next(x0 = x_ini_start_mhe, 
                     p = vertcat(par_val, opt_u_nom[:,i]),
                     z0 = z_ini_start_mhe)

    J_cl = J_cl  + u7_coef*opt_u_nom[7,i] + \
               (u8_coef2[i]*opt_u_nom[8,i] + u9_coef*(opt_u_nom[9,i]+opt_u_nom[10,i]+opt_u_nom[11,i]))

    J_cl_arc_power = J_cl_arc_power + u8_coef2[i]*opt_u_nom[8,i]

    # Electricity price change 
    if i == n_change-1:
         u8_coef1[n_change:] = u8_coef1_updated
         u8_coef = [x1_ * (1/60) for x1_ in u8_coef1]          
    #if i == n_mhe-1: 
    #    J_cl = J_cl - xk15o_coef*((Ik_next['xf'][15]*xk15i_coef*xk15g_coef) + (((Ik_next['zf'][42])*exp(Ik_next['zf'][44]))*xk15i_coef))
                      
    # Fix these for the next measurment calculations    
    x_ini_start_mhe = Ik_next['xf']
    z_ini_start_mhe = Ik_next['zf'] 
    
    # Simulate the model 

    # Build the state matrices
    x_mat_act = horzcat(x_mat_act, Ik_next['xf'])
    z_mat_act = horzcat(z_mat_act, Ik_next['zf'])
    
    # Get all the dependent variable
    d_all_next = f_y(Ik_next['xf'], Ik_next['zf'], opt_u_nom[:,i], 
                     par_val) 
   
    # Choose to get the measurement vector
    y_next = vertcat(d_all_next[318:321], d_all_next[322], 
                     Ik_next['xf'][26:28],               # continuous measurements
                     Ik_next['xf'][19], d_all_next[396], # molten metal measurements 
                     d_all_next[196], d_all_next[201], d_all_next[203], 
                     d_all_next[205:207])                # slag measurements          
    
    # Build the measurement matrix
    y_mat = horzcat(y_mat, y_next)
    
    # Add noise to the measurements
    NP.random.seed(i+1)
    y_next_wn = y_next + sigma_y_full*NP.random.randn(y_next.shape[0],1) 

    # Build the measurement matrix (with noise)
    y_mat_wn = horzcat(y_mat_wn, y_next_wn)    
    
    #%% NLP formulation    
    # Start with an empty NLP
    w=[]
    w0 = []
    lbw = []
    ubw = []
    g=[]
    lbg = []
    ubg = []    
    mhe_obj = 0

    # Declare the end measurement parameter
    Nmeas_end = meas_struct[i+1]
    p_nlp0 = y_mat_wn[0:Nmeas_end, i+1]

    # "Lift" initial conditions
    X_ini_def = SX.sym('X_ini_def', Nstates)
    w += [X_ini_def]
    w0 += X_ini.nonzeros() # noise added so take care

    if disturb_add ==1:
        for i5 in range(Nstates):
            if i5 in [0,1,2,3,4,5,6,7,8,9,10,11,12 , 14,15,16,17,18,19,20,21,22,23,24,25,26,27]:#range(Nstates-1):#[26,27,7]:#,7,8,10]:
                if abs(X_ini[i5]) > 0 and abs(X_ini[i5]) < 5:
                    if X_ini[i5] < 0:
                        lbw += [-10]
                    else:
                        lbw += [0]
                else:
                    if X_ini[i5] <= 0:
                        lbw += (1.5*X_ini[i5]).nonzeros()
                    else:    
                        lbw += (0.5*X_ini[i5]).nonzeros()
            else:
                if i5 ==13:
                    lbw += (0.8*X_ini[i5]).nonzeros()
                else:
                    lbw += [-0.1]
                #lbw += [X_ini.nonzeros()[i5]] 
        
        for i4 in range(Nstates):
            if i4 in [0,1,2,3,4,5,6,7,8,9,10,11,12 , 14,15,16,17,18,19,20,21,22,23,24,25,26,27]:#range(Nstates-1):#[26,27,7]:#,7,8,10]:          
                if abs(X_ini[i4]) > 0 and abs(X_ini[i4]) < 5:
                    if X_ini[i4] < 0:
                        ubw += [0]
                    else:
                        ubw += [10]
                else:
                    if X_ini[i4] < 0:
                        if X_ini[i4] > -0.5: # unwanted negative values due to noise additions
                            ubw += 0.5 # fixed upper bound 
                        else:    
                            ubw += (0.5*X_ini[i4]).nonzeros()
                    else:    
                        #ubw += (1.5*X_ini[i4]+0.0001).nonzeros()
                        ubw += (1.5*X_ini[i4]+0).nonzeros()
            else:
                if i4 == 13:
                    ubw +=(1.2*X_ini[i4]).nonzeros()
                else:    
                    ubw += [0.1]
                #ubw += [X_ini.nonzeros()[i4]] 
    else:          
        for i5 in range(Nstates):
            if i5 in range(Nstates-0):#[26,27,7]:#,7,8,10]:
                if abs(X_ini[i5]) > 0 and abs(X_ini[i5]) < 5:
                    if X_ini[i5] < 0:
                        lbw += [-10]
                    else:
                        lbw += [0]
                else:
                    if X_ini[i5] <= 0:
                        lbw += (1.5*X_ini[i5]).nonzeros()
                    else:    
                        lbw += (0.5*X_ini[i5]).nonzeros()
            else:
                #lbw += [-0.1]
                lbw += [X_ini.nonzeros()[i5]] 
        
        for i4 in range(Nstates):
            if i4 in range(Nstates-0):#[26,27,7]:#,7,8,10]:            
                if abs(X_ini[i4]) > 0 and abs(X_ini[i4]) < 5:
                    if X_ini[i4] < 0:
                        ubw += [0]
                    else:
                        ubw += [10]
                else:
                    if X_ini[i4] < 0:
                        if X_ini[i4] > -0.5: # unwanted negative values due to noise additions
                            ubw += 0.5 # fixed upper bound 
                        else:    
                            ubw += (0.5*X_ini[i4]).nonzeros()
                    else:    
                        ubw += (1.5*X_ini[i4]+0.0001).nonzeros()  
            else:
                #ubw += [0.1]
                ubw += [X_ini.nonzeros()[i4]] 
                
    #"""        
    # Arrival cost term of mhe obj function
    mhe_obj += mtimes([(X_ini_def - X_ini).T, S, (X_ini_def - X_ini)])
    
    # Add the initial algebraic variables as decision variables 
    Z_ini_def = SX.sym('Z_ini_def', Nalgvars)
    w += [Z_ini_def]
    w0 += Z_ini.nonzeros()

    # Bounds for initial algebraic variables        
    for i3 in range(Nalgvars):
        if i3 in range(14)+[15,17,18,20,22,24,26]+range(28,43)+range(52,108)+[110,113,115,119,120]: 
            if abs(Z_ini[i3]) > 0 and abs(Z_ini[i3]) < 5:
                if Z_ini[i3] < 0:
                    lbw += [-10]
                else:
                    lbw += [0]
            else:
                if Z_ini[i3] <= 0:
                    lbw += (1.5*Z_ini[i3]).nonzeros()
                else:    
                    lbw += (0.5*Z_ini[i3]).nonzeros()
        #else:      #  0     1   2   3    4    5   6   7    8   9   10  11   12  13   14   15     
        #    if i3 in [14*, 51, 16, 44*, 43*, 19, 45, 21*, 46, 23*, 47, 25*, 48, 27*, 49*, 50*]
        else:      #  0     1   2   3    4    5   6   7    8   9   10  11   12  13   14   15     
            if i3 in [51, 16, 19, 45, 46, 47, 48]:
                lbw += (1.2*Z_ini[i3]).nonzeros()
            else:
                lbw += (1.5*Z_ini[i3]).nonzeros()
                #lbw += [-1e10]
        
    for i2 in range(Nalgvars):
        if i2 in range(14)+[15,17,18,20,22,24,26]+range(28,43)+range(52,108)+[110,113,115,119,120]: 
            if abs(Z_ini[i2]) > 0 and abs(Z_ini[i2]) < 5:
                if Z_ini[i2] < 0:
                    ubw += [0]
                else:
                    ubw += [10]
            else:
                if Z_ini[i2] < 0:
                    ubw += (0.5*Z_ini[i2]).nonzeros()
                else:    
                    ubw += (1.5*Z_ini[i2]+0.0001).nonzeros()
        else:      #  0     1   2   3    4    5   6   7    8   9   10  11   12  13   14   15     
            if i2 in [51, 16, 19, 45, 46, 47, 48]:
                ubw += (0.8*Z_ini[i2]).nonzeros()
            else:
                #ubw += (0.9*Z_ini_be[i2]).nonzeros()                    
                #ubw += [0.001] # [0.0001] upto 47      
                ubw += [0.0001]# upto 47 
            
    # Backward euler related equation                   
    Xk = X_ini
    
    #######################################################    
    ## Model and measurement noise noise addition loop
    #######################################################
    
    #%% Condition to check if batch MHE or not
    if i >= mhe_horizon: #(non-batch MHE)

        # Get the measurement structure 
        Nmeas = meas_struct[i-mhe_horizon+1:i+2]
        
        for e1 in range(mhe_horizon):
            
            if e1 == mhe_horizon-1: # if last element then 
                # Initialization of backeuler simulation
                # Try this
                Ik1 = I_full(x0=X_ini_be, 
                             p=vertcat(par_val_mismatch, 
                                   opt_u_nom[:,i+n_ini]),
                                   z0 = Z_ini_be)
                  
            else:
                # Use from previous solve
                x_guess1 = x_guess[:,steps*e1:steps+steps*e1]
                z_guess1 = z_guess[:,steps*e1:steps+steps*e1]            
            
            # New NLP variable for the controls (Inputs + model noises)
            Uk = SX.sym('U_' + str(i) + '_' + str(e1), Ninputs + Nstates)
            w   += [Uk]
    
            # Fix inputs
            lbw += [k for k in opt_u_nom[0:Ninputs,i - (mhe_horizon -1) + e1 + n_ini].nonzeros()]
            ubw += [k for k in opt_u_nom[0:Ninputs,i - (mhe_horizon -1) + e1 + n_ini].nonzeros()]
            
            #lbw += U_dm2[i - (mhe_horizon -1) + e1 + n_ini] + U_dm1[i - (mhe_horizon -1) + e1 + n_ini]
            #ubw += U_dm2[i - (mhe_horizon -1) + e1 + n_ini] + U_dm1[i - (mhe_horizon -1) + e1 + n_ini]

            # Provide bounds for model noises 
            lbw += [k for k in w_lb.nonzeros()]
            ubw += [k for k in w_ub.nonzeros()]
            #lbw += [-1]*Nstates    
            #ubw += [1]*Nstates

            # Initial guess for model noises
            w0 += [k for k in opt_u_nom[0:Ninputs,i - (mhe_horizon -1) + e1 + n_ini].nonzeros()]
            
            #w0 += U_dm2[i - (mhe_horizon -1) + e1 + n_ini] + U_dm1[i - (mhe_horizon -1) + e1 + n_ini]
            if e1 == mhe_horizon-1: # if last element
                #w0 += [0]*Nstates 
                # Instead of using zeros, we initialize w  with the previous 
                # time step value
                w0 += [k for k in u_opt[Ninputs:Ninputs+Nstates,e1-1].nonzeros()]
            else:     
                # pass the solution from last solve
                w0 += [k for k in u_opt[Ninputs:Ninputs+Nstates,e1].nonzeros()]                
        
            # Add contribution of model noises to cost function
            mhe_obj += mtimes([(Uk[12:]).T, Q, (Uk[12:])])
            
            #%% Add initial measurement noise to mhe objective
            if e1 == 0:
                # Add algebraic equations at t=0 as equality constraints
                fj, zj = f(X_ini_def, Z_ini_def, par_val_mismatch, 
                           vertcat(opt_u_nom[0:Ninputs,i-mhe_horizon+1+n_ini], Uk[12:]))       
                g += [zj]
                lbg += [0]*Nalgvars
                ubg += [0]*Nalgvars
                
                # First we have to get dependent vars from alg vars
                d_first = f_y(X_ini_def, Z_ini_def, 
                              vertcat(opt_u_nom[0:Ninputs,i-mhe_horizon+1+n_ini], Uk[12:]), 
                              par_val_mismatch) 
                
                # Choose your measurements
                y_0 = vertcat(d_first[318:321], d_first[322], 
                              X_ini_def[26:28],               # continuous measurements
                              X_ini_def[19], d_first[396], # molten metal measurements 
                              d_first[196], d_first[201], d_first[203], 
                              d_first[205:207])                #slag measurements 
                 
                # Now add the cost to mhe objective function 
                if Nmeas[0] == 13:
                    mhe_obj += mtimes([(y_0[0:Nmeas[0]]-y_mat_wn[0:Nmeas[0],i-mhe_horizon+1]).T,
                                       R_full,
                                       (y_0[0:Nmeas[0]]-y_mat_wn[0:Nmeas[0],i-mhe_horizon+1])])
                                   
                if Nmeas[0] == 8:
                    mhe_obj += mtimes([(y_0[0:Nmeas[0]]-y_mat_wn[0:Nmeas[0],i-mhe_horizon+1]).T,
                                       R_half,
                                       (y_0[0:Nmeas[0]]-y_mat_wn[0:Nmeas[0],i-mhe_horizon+1])])

                if Nmeas[0] == 6:
                    mhe_obj += mtimes([(y_0[0:Nmeas[0]]-y_mat_wn[0:Nmeas[0],i-mhe_horizon+1]).T,
                                       R,
                                       (y_0[0:Nmeas[0]]-y_mat_wn[0:Nmeas[0],i-mhe_horizon+1])])              
            
            #%% Further Measurement noise additions (all except for t=0)
            for j in range(steps):

                # Initialization of backward euler
                if e1 == mhe_horizon-1: 
                    X_ini_be = Ik1['xf'][:,j+1]
                    Z_ini_be = Ik1['zf'][:,j+1]
                else:
                    X_ini_be = x_guess1[:,j]
                    Z_ini_be = z_guess1[:,j]
                
                # State and alg vars at the euler points
                Xk_end = SX.sym('X_' + str(i) + '_' + str(e1+1) + '_' + str(j),
                                Nstates)
                Zk_end = SX.sym('Zz_' + str(i) + '_' + str(e1+1) + '_' + str(j),
                                Nalgvars)
                
                # Both states and alg vars are decision variables
                w += [Xk_end]
                w += [Zk_end]        

                # Initial guesses
                w0 += [k for k in X_ini_be.nonzeros()]
                w0 += [k for k in Z_ini_be.nonzeros()]

                if disturb_add ==1:
                    for i5 in range(Nstates):
                        if i5 in [0,1,2,3,4,5,6,7,8,9,10,11,12 , 14,15,16,17,18,19,20,21,22,23,24,25,26,27]:#range(Nstates-1):#[26,27,7]:#,7,8,10]:
                            if abs(X_ini_be[i5]) > 0 and abs(X_ini_be[i5]) < 5:
                                if X_ini_be[i5] < 0:
                                    lbw += [-10]
                                else:
                                    lbw += [0]
                            else:
                                if X_ini_be[i5] <= 0:
                                    lbw += (1.5*X_ini_be[i5]).nonzeros()
                                else:    
                                    lbw += (0.5*X_ini_be[i5]).nonzeros()
                        else:
                            if i5 ==13:
                                lbw += (0.8*X_ini_be[i5]).nonzeros()
                            else:
                                lbw += [-0.1]
                            #lbw += [X_ini_be.nonzeros()[i5]] 
                    
                    for i4 in range(Nstates):
                        if i4 in [0,1,2,3,4,5,6,7,8,9,10,11,12 , 14,15,16,17,18,19,20,21,22,23,24,25,26,27]:#range(Nstates-1):#[26,27,7]:#,7,8,10]:          
                            if abs(X_ini_be[i4]) > 0 and abs(X_ini_be[i4]) < 5:
                                if X_ini_be[i4] < 0:
                                    ubw += [0]
                                else:
                                    ubw += [10]
                            else:
                                if X_ini_be[i4] < 0:
                                    if X_ini_be[i4] > -0.5: # unwanted negative values due to noise additions
                                        ubw += 0.5 # fixed upper bound 
                                    else:    
                                        ubw += (0.5*X_ini_be[i4]).nonzeros()
                                else:    
                                    #ubw += (1.5*X_ini_be[i4]+0.0001).nonzeros()
                                    ubw += (1.5*X_ini_be[i4]+0).nonzeros()
                        else:
                            if i4 == 13:
                                ubw +=(1.5*X_ini_be[i4]).nonzeros()
                            else:    
                                ubw += [0.1]
                            #ubw += [X_ini_be.nonzeros()[i4]] 
                else:          
                    for i5 in range(Nstates):
                        if i5 in range(Nstates-0):#[26,27,7]:#,7,8,10]:
                            if abs(X_ini_be[i5]) > 0 and abs(X_ini_be[i5]) < 5:
                                if X_ini_be[i5] < 0:
                                    lbw += [-10]
                                else:
                                    lbw += [0]
                            else:
                                if X_ini_be[i5] <= 0:
                                    lbw += (1.5*X_ini_be[i5]).nonzeros()
                                else:    
                                    lbw += (0.5*X_ini_be[i5]).nonzeros()
                        else:
                            #lbw += [-0.1]
                            lbw += [X_ini_be.nonzeros()[i5]] 
                    
                    for i4 in range(Nstates):
                        if i4 in range(Nstates-0):#[26,27,7]:#,7,8,10]:            
                            if abs(X_ini_be[i4]) > 0 and abs(X_ini_be[i4]) < 5:
                                if X_ini_be[i4] < 0:
                                    ubw += [0]
                                else:
                                    ubw += [10]
                            else:
                                if X_ini_be[i4] < 0:
                                    if X_ini_be[i4] > -0.5: # unwanted negative values due to noise additions
                                        ubw += 0.5 # fixed upper bound 
                                    else:    
                                        ubw += (0.5*X_ini_be[i4]).nonzeros()
                                else:    
                                    ubw += (1.5*X_ini_be[i4]+0.0001).nonzeros()  
                        else:
                            #ubw += [0.1]
                            ubw += [X_ini_be.nonzeros()[i4]] 

                # Bounds for algebraic variables        
                for i3 in range(Nalgvars):
                    if i3 in range(14)+[15,17,18,20,22,24,26]+range(28,43)+range(52,108)+[110,113,115,119,120]: 
                        if abs(Z_ini_be[i3]) > 0 and abs(Z_ini_be[i3]) < 5:
                            if Z_ini_be[i3] < 0:
                                lbw += [-10]
                            else:
                                lbw += [0]
                        else:
                            if Z_ini_be[i3] <= 0:
                                lbw += (1.5*Z_ini_be[i3]).nonzeros()
                            else:    
                                lbw += (0.5*Z_ini_be[i3]).nonzeros()
                    #else:      #  0     1   2   3    4    5   6   7    8   9   10  11   12  13   14   15     
                    #    if i3 in [14*, 51, 16, 44*, 43*, 19, 45, 21*, 46, 23*, 47, 25*, 48, 27*, 49*, 50*]
                    else:      #  0     1   2   3    4    5   6   7    8   9   10  11   12  13   14   15     
                        if i3 in [51, 16, 19, 45, 46, 47, 48]:
                            lbw += (1.2*Z_ini_be[i3]).nonzeros()
                        else:
                            lbw += (1.5*Z_ini_be[i3]).nonzeros()
                            #lbw += [-1e10]
        
                for i2 in range(Nalgvars):
                    if i2 in range(14)+[15,17,18,20,22,24,26]+range(28,43)+range(52,108)+[110,113,115,119,120]: 
                        if abs(Z_ini_be[i2]) > 0 and abs(Z_ini_be[i2]) < 5:
                            if Z_ini_be[i2] < 0:
                                ubw += [0]
                            else:
                                ubw += [10]
                        else:
                            if Z_ini_be[i2] < 0:
                                ubw += (0.5*Z_ini_be[i2]).nonzeros()
                            else:    
                                ubw += (1.5*Z_ini_be[i2]+0.0001).nonzeros()
                    else:      #  0     1   2   3    4    5   6   7    8   9   10  11   12  13   14   15     
                        if i2 in [51, 16, 19, 45, 46, 47, 48]:
                            ubw += (0.8*Z_ini_be[i2]).nonzeros()
                        else:
                            #ubw += (0.9*Z_ini_be[i2]).nonzeros()                    
                            #ubw += [0.001] # [0.0001] upto 47      
                            ubw += [0.0001]# upto 47                       
                #lbw += [k for k in Z_ini_be.nonzeros()]
                #ubw += [k for k in Z_ini_be.nonzeros()]  
            
                # Euler equations
                fj, zj = f(Xk_end, Zk_end, par_val_mismatch, Uk)       
                g += [Xk - Xk_end + fj*h]
                lbg += [0]*Nstates
                ubg += [0]*Nstates
                
                g += [zj]
                lbg += [0]*Nalgvars
                ubg += [0]*Nalgvars
            
                # Measurement noise contribution (at the end of the finite element)
                if j == steps-1:  # If last step in the finite element

                    # First we have to get dependent vars from alg vars
                    d_mid = f_y(Xk_end, Zk_end, Uk, par_val_mismatch)

                    # Choose your measurements
                    y_mid = vertcat(d_mid[318:321], d_mid[322],
                                    Xk_end[26:28],               # continuous measurements
                                    Xk_end[19], d_mid[396],      # molten metal measurements
                                    d_mid[196], d_mid[201], d_mid[203],
                                    d_mid[205:207])                # slag measurements

                    # Now add the cost to mhe objective function
                    if e1 == mhe_horizon-1:  # If last element use the measurement parameter
                        if Nmeas[e1 + 1] == 13:
                            mhe_obj += mtimes([(p_nlp0 - y_mat_wn[0:Nmeas[e1 + 1],
                                                                         i - mhe_horizon + 1 + e1 + 1]).T,
                                               R_full,
                                               (p_nlp0 - y_mat_wn[0:Nmeas[e1 + 1],
                                                                         i - mhe_horizon + 1 + e1 + 1])])

                        if Nmeas[e1 + 1] == 8:
                            mhe_obj += mtimes([(p_nlp0 - y_mat_wn[0:Nmeas[e1 + 1],
                                                                         i - mhe_horizon + 1 + e1 + 1]).T,
                                               R_half,
                                               (p_nlp0 - y_mat_wn[0:Nmeas[e1 + 1],
                                                                         i - mhe_horizon + 1 + e1 + 1])])

                        if Nmeas[e1 + 1] == 6:
                            mhe_obj += mtimes([(p_nlp0 - y_mat_wn[0:Nmeas[e1 + 1],
                                                                         i - mhe_horizon + 1 + e1 + 1]).T,
                                               R,
                                               (p_nlp0 - y_mat_wn[0:Nmeas[e1 + 1],
                                                                         i - mhe_horizon + 1 + e1 + 1])])
                    else:
                        if Nmeas[e1+1] == 13:
                            mhe_obj += mtimes([(y_mid[0:Nmeas[e1+1]]-y_mat_wn[0:Nmeas[e1+1],
                                                                     i-mhe_horizon+1+e1+1]).T,
                                               R_full,
                                               (y_mid[0:Nmeas[e1+1]]-y_mat_wn[0:Nmeas[e1+1],
                                                                     i-mhe_horizon+1+e1+1])])

                        if Nmeas[e1+1] == 8:
                            mhe_obj += mtimes([(y_mid[0:Nmeas[e1+1]]-y_mat_wn[0:Nmeas[e1+1],
                                                                     i-mhe_horizon+1+e1+1]).T,
                                               R_half,
                                               (y_mid[0:Nmeas[e1+1]]-y_mat_wn[0:Nmeas[e1+1],
                                                                     i-mhe_horizon+1+e1+1])])

                        if Nmeas[e1+1] == 6:
                            mhe_obj += mtimes([(y_mid[0:Nmeas[e1+1]]-y_mat_wn[0:Nmeas[e1+1],
                                                                     i-mhe_horizon+1+e1+1]).T,
                                               R,
                                               (y_mid[0:Nmeas[e1+1]]-y_mat_wn[0:Nmeas[e1+1],
                                                                     i-mhe_horizon+1+e1+1])])

                # Add equality constraint
                Xk = Xk_end
    
    #%% if batch MHE                               
    else: 
        # Get the measurement structure for the horizon
        Nmeas = meas_struct[0:i+2]        
        
        for e1 in range(i+1):
                        
            if e1 == i: # if last element then 
                # Initialization of backeuler simulation
                Ik1 = I_full(x0=X_ini_be, 
                                 p=vertcat(par_val_mismatch, opt_u_nom[:,e1+n_ini]), 
                                           z0 = Z_ini_be)                                 
            else:
                # Use from previous solve
                x_guess1 = x_guess[:,steps*e1:steps+steps*e1]
                z_guess1 = z_guess[:,steps*e1:steps+steps*e1]
            
            # New NLP variable for the controls (Inputs + model noises)
            Uk = SX.sym('U_' + str(i) + '_' + str(e1), Ninputs + Nstates)
            w   += [Uk]
    
            # Fix inputs 
            lbw += [k for k in opt_u_nom[0:Ninputs, e1 + n_ini].nonzeros()]
            ubw += [k for k in opt_u_nom[0:Ninputs, e1 + n_ini].nonzeros()]

            # Provide bounds for model noises
            lbw += [k for k in w_lb.nonzeros()]
            ubw += [k for k in w_ub.nonzeros()]   
            #lbw += [-1]*Nstates    
            #ubw += [1]*Nstates
        
            # Initial guess
            w0 += [k for k in opt_u_nom[0:Ninputs, e1 + n_ini].nonzeros()]
            if e1 == i: # if last element
                w0 += [0]*Nstates 
            else:     
                # pass the solution from last solve
                if i == 0: # very first solve is a vector not a matrix
                    w0 += [k for k in u_opt[Ninputs:Ninputs+Nstates].nonzeros()]
                else:    
                    w0 += [k for k in u_opt[Ninputs:Ninputs+Nstates,e1].nonzeros()]
        
            # Add contribution of model noises to cost function
            mhe_obj += mtimes([(Uk[12:]).T, Q, (Uk[12:])])
            
            #%% Add initial measurement noise to mhe objective
            if e1 == 0:
                # Add algebraic equations at t=0 as equality constraints
                fj, zj = f(X_ini_def, Z_ini_def, par_val_mismatch, 
                           vertcat(opt_u_nom[0:Ninputs, n_ini], Uk[12:])) 
                           
                g += [zj]
                lbg += [0]*Nalgvars
                ubg += [0]*Nalgvars
                
                # First we have to get dependent vars from alg vars
                d_first = f_y(X_ini_def, Z_ini_def, 
                              vertcat(opt_u_nom[0:Ninputs, n_ini], Uk[12:]), 
                              par_val_mismatch) 
                
                # Choose your measurements
                y_0 = vertcat(d_first[318:321], d_first[322], 
                              X_ini_def[26:28],               # continuous measurements
                              X_ini_def[19], d_first[396], # molten metal measurements 
                              d_first[196], d_first[201], d_first[203], 
                              d_first[205:207])                #slag measurements 
                 
                # Now add the cost to mhe objective function
                 
                if Nmeas[0] == 13:
                    mhe_obj += mtimes([(y_0[0:Nmeas[0]]-y_mat_wn[0:Nmeas[0],0]).T,
                                       R_full,
                                       (y_0[0:Nmeas[0]]-y_mat_wn[0:Nmeas[0],0])])
                                   
                if Nmeas[0] == 8:
                    mhe_obj += mtimes([(y_0[0:Nmeas[0]]-y_mat_wn[0:Nmeas[0],0]).T,
                                       R_half,
                                       (y_0[0:Nmeas[0]]-y_mat_wn[0:Nmeas[0],0])])

                if Nmeas[0] == 6:
                    mhe_obj += mtimes([(y_0[0:Nmeas[0]]-y_mat_wn[0:Nmeas[0],0]).T,
                                       R,
                                       (y_0[0:Nmeas[0]]-y_mat_wn[0:Nmeas[0],0])])                      
                
            #%% Further Measurement noise additions (all except for t=0)
            for j in range(steps):

                # Initialization of backward euler
                if e1 == i: 
                    X_ini_be = Ik1['xf'][:,j+1]
                    Z_ini_be = Ik1['zf'][:,j+1]
                else:
                    X_ini_be = x_guess1[:,j]
                    Z_ini_be = z_guess1[:,j]
                
                # State and alg vars at the euler points
                Xk_end = SX.sym('X_' + str(i) + '_' + str(e1+1) + '_' + str(j),
                                Nstates)
                Zk_end = SX.sym('Zz_' + str(i) + '_' + str(e1+1) + '_' + str(j),
                                Nalgvars)
                
                # Both states and alg vars are decision variables
                w += [Xk_end]
                w += [Zk_end]        

                # Initial guesses
                w0 += [k for k in X_ini_be.nonzeros()]
                w0 += [k for k in Z_ini_be.nonzeros()]
                
                # Bounds for state variables
                if disturb_add ==1:
                    for i5 in range(Nstates):
                        if i5 in range(Nstates-1):#[26,27,7]:#,7,8,10]:
                            if abs(X_ini_be[i5]) > 0 and abs(X_ini_be[i5]) < 5:
                                if X_ini_be[i5] < 0:
                                    lbw += [-10]
                                else:
                                    lbw += [0]
                            else:
                                if X_ini_be[i5] <= 0:
                                    lbw += (1.5*X_ini_be[i5]).nonzeros()
                                else:    
                                    lbw += (0.5*X_ini_be[i5]).nonzeros()
                        else: 
                            lbw += [-0.1]
                            #lbw+= [X_ini_be.nonzeros()[i5]] 
                       
                    for i4 in range(Nstates):
                        if i4 in range(Nstates-1):#[26,27,7]:#,7,8,10]:            
                            if abs(X_ini_be[i4]) > 0 and abs(X_ini_be[i4]) < 5:
                                if X_ini_be[i4] < 0:
                                    ubw += [0]
                                else:
                                    ubw += [10]
                            else:
                                if X_ini_be[i4] < 0:
                                    ubw += (0.5*X_ini_be[i4]).nonzeros()
                                else:    
                                    ubw += (1.5*X_ini_be[i4]+0.0001).nonzeros()  
                        else:
                            ubw += [0.1]
                            #ubw+= [X_ini_be.nonzeros()[i4]] 
                else:                                
                    for i5 in range(Nstates):
                        if i5 in range(Nstates-0):#[26,27,7]:#,7,8,10]:
                            if abs(X_ini_be[i5]) > 0 and abs(X_ini_be[i5]) < 5:
                                if X_ini_be[i5] < 0:
                                    lbw += [-10]
                                else:
                                    lbw += [0]
                            else:
                                if X_ini_be[i5] <= 0:
                                    lbw += (1.5*X_ini_be[i5]).nonzeros()
                                else:    
                                    lbw += (0.5*X_ini_be[i5]).nonzeros()
                        else: 
                            #lbw += [-0.1]
                            lbw+= [X_ini_be.nonzeros()[i5]] 
                       
                    for i4 in range(Nstates):
                        if i4 in range(Nstates-0):#[26,27,7]:#,7,8,10]:            
                            if abs(X_ini_be[i4]) > 0 and abs(X_ini_be[i4]) < 5:
                                if X_ini_be[i4] < 0:
                                    ubw += [0]
                                else:
                                    ubw += [10]
                            else:
                                if X_ini_be[i4] < 0:
                                    ubw += (0.5*X_ini_be[i4]).nonzeros()
                                else:    
                                    ubw += (1.5*X_ini_be[i4]+0.0001).nonzeros()  
                        else:
                            #ubw += [0.1]
                            ubw+= [X_ini_be.nonzeros()[i4]] 
                                                        
                # Bounds for algebraic variables        
                for i3 in range(Nalgvars):
                    if i3 in range(14)+[15,17,18,20,22,24,26]+range(28,43)+range(52,108)+[110,113,115,119,120]: 
                        if abs(Z_ini_be[i3]) > 0 and abs(Z_ini_be[i3]) < 5:
                            if Z_ini_be[i3] < 0:
                                lbw += [-10]
                            else:
                                lbw += [0]
                        else:
                            if Z_ini_be[i3] <= 0:
                                lbw += (1.5*Z_ini_be[i3]).nonzeros()
                            else:    
                                lbw += (0.5*Z_ini_be[i3]).nonzeros()
                    #else:      #  0     1   2   3    4    5   6   7    8   9   10  11   12  13   14   15     
                    #    if i3 in [14*, 51, 16, 44*, 43*, 19, 45, 21*, 46, 23*, 47, 25*, 48, 27*, 49*, 50*]
                    else:      #  0     1   2   3    4    5   6   7    8   9   10  11   12  13   14   15     
                        if i3 in [51, 16, 19, 45, 46, 47, 48]:
                            lbw += (1.2*Z_ini_be[i3]).nonzeros()
                        else:
                            lbw += (1.5*Z_ini_be[i3]).nonzeros()
                            #lbw += [-1e10]
        
                for i2 in range(Nalgvars):
                    if i2 in range(14)+[15,17,18,20,22,24,26]+range(28,43)+range(52,108)+[110,113,115,119,120]: 
                        if abs(Z_ini_be[i2]) > 0 and abs(Z_ini_be[i2]) < 5:
                            if Z_ini_be[i2] < 0:
                                ubw += [0]
                            else:
                                ubw += [10]
                        else:
                            if Z_ini_be[i2] < 0:
                                ubw += (0.5*Z_ini_be[i2]).nonzeros()
                            else:    
                                ubw += (1.5*Z_ini_be[i2]+0.0001).nonzeros()
                    else:      #  0     1   2   3    4    5   6   7    8   9   10  11   12  13   14   15     
                        if i2 in [51, 16, 19, 45, 46, 47, 48]:
                            ubw += (0.8*Z_ini_be[i2]).nonzeros()
                        else:
                            #ubw += (0.9*Z_ini_be[i2]).nonzeros()                    
                            #ubw += [0.001] # [0.0001] upto 47      
                            ubw += [0.0001]# upto 47                       
                #lbw += [k for k in Z_ini_be.nonzeros()]
                #ubw += [k for k in Z_ini_be.nonzeros()]  
            
                # Euler equations
                fj, zj = f(Xk_end, Zk_end, par_val_mismatch, Uk)       
                g += [Xk - Xk_end + fj*h]
                lbg += [0]*Nstates
                ubg += [0]*Nstates
                
                g += [zj]
                lbg += [0]*Nalgvars
                ubg += [0]*Nalgvars
            
                if j == steps-1:  # If last step in the finite element

                    # First we have to get dependent vars from alg vars
                    d_mid = f_y(Xk_end, Zk_end, Uk, par_val_mismatch)

                    # Choose your measurements
                    y_mid = vertcat(d_mid[318:321], d_mid[322],
                                    Xk_end[26:28],               # continuous measurements
                                    Xk_end[19], d_mid[396],      # molten metal measurements
                                    d_mid[196], d_mid[201], d_mid[203],
                                    d_mid[205:207])                # slag measurements

                    # Now add the cost to mhe objective function
                    if e1 == i:  # if last element then use measurement parameter declaration
                        if Nmeas[e1 + 1] == 13:
                            mhe_obj += mtimes([(p_nlp0 - y_mat_wn[0:Nmeas[e1 + 1],
                                                                         e1 + 1]).T,
                                               R_full,
                                               (p_nlp0 - y_mat_wn[0:Nmeas[e1 + 1],
                                                                         e1 + 1])])

                        if Nmeas[e1 + 1] == 8:
                            mhe_obj += mtimes([(p_nlp0 - y_mat_wn[0:Nmeas[e1 + 1],
                                                                         e1 + 1]).T,
                                               R_half,
                                               (p_nlp0 - y_mat_wn[0:Nmeas[e1 + 1],
                                                                         e1 + 1])])

                        if Nmeas[e1 + 1] == 6:
                            mhe_obj += mtimes([(p_nlp0 - y_mat_wn[0:Nmeas[e1 + 1],
                                                                         e1 + 1]).T,
                                               R,
                                               (p_nlp0 - y_mat_wn[0:Nmeas[e1 + 1],
                                                                         e1 + 1])])
                    else:
                        if Nmeas[e1+1] == 13:
                            mhe_obj += mtimes([(y_mid[0:Nmeas[e1+1]]-y_mat_wn[0:Nmeas[e1+1],
                                                                     e1+1]).T,
                                               R_full,
                                               (y_mid[0:Nmeas[e1+1]]-y_mat_wn[0:Nmeas[e1+1],
                                                                     e1+1])])

                        if Nmeas[e1+1] == 8:
                            mhe_obj += mtimes([(y_mid[0:Nmeas[e1+1]]-y_mat_wn[0:Nmeas[e1+1],
                                                                     e1+1]).T,
                                               R_half,
                                               (y_mid[0:Nmeas[e1+1]]-y_mat_wn[0:Nmeas[e1+1],
                                                                     e1+1])])

                        if Nmeas[e1+1] == 6:
                            mhe_obj += mtimes([(y_mid[0:Nmeas[e1+1]]-y_mat_wn[0:Nmeas[e1+1],
                                                                     e1+1]).T,
                                               R,
                                               (y_mid[0:Nmeas[e1+1]]-y_mat_wn[0:Nmeas[e1+1],
                                                                     e1+1])])

                # Add equality constraint
                Xk = Xk_end                                 
    
    #%%#################################################
    #            Solve the MHE problem
    ####################################################
    
    # Print statements
    print '\n****************************************************************'
    print '                       Iteration No. '+ repr(i)
    print '****************************************************************\n'
    # NLP solver options
    opts_nlp = {}
    opts_nlp["ipopt.linear_solver"] = 'ma27'
    opts_nlp["ipopt.max_iter"] = 1000

    if i > 0: # Ue warmstart when past initial iteration
        opts_nlp["ipopt.warm_start_init_point"] = 'yes'
        opts_nlp["ipopt.warm_start_bound_push"] = 1e-6
        opts_nlp["ipopt.warm_start_mult_bound_push"] = 1e-6
        opts_nlp["ipopt.mu_init"] = 1e-6    
    #opts_nlp["ipopt.linear_system_scaling"] = 'slack-based'
    #opts_nlp["ipopt.alpha_for_y"]='primal-and-full'
    #opts_nlp["ipopt.obj_scaling_factor"] = 1e-3
    #opts_nlp["ipopt.hessian_approximation"] = 'limited-memory'
    prob = {'f': mhe_obj, 'x': vertcat(*w), 'g': vertcat(*g)}
    solver = nlpsol('solver', 'ipopt', prob, opts_nlp)
        
    # Warm start the nlp solve with dual solutions from previous solve
    zero_padding = DM.zeros(Ninputs + Nstates*(steps+1) + Nalgvars*(steps),1)
    zero_padding1 = DM.zeros(Nstates*(steps) + Nalgvars*(steps))
				
    # Start clock to measure simulation time
    startTime = time.time()
				
    if i <= mhe_horizon-1:
        if i == 0:  # if first solve no dual guess available
            sol1 = solver(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg)#, p=p_nlp0)
        else:  # pass dual and add zeros in the end
            guess2 = vertcat(w_opt_dual, zero_padding)
            guess3 = vertcat(w_opt_lam_g, zero_padding1)

            sol1 = solver(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg,
                         lam_x0=guess2, lam_g0=guess3) #, p=p_nlp0)
    else:  # nonbatch mhe
        # Drop start and add zeros in the end
        guess2 = vertcat(w_opt_dual[Nstates*(steps+1) + Ninputs + Nalgvars * steps:],
                         zero_padding)
        # Drop start and add zeros in the end
        guess3 = vertcat(w_opt_lam_g[Nstates * steps + Nalgvars * steps:],
                         zero_padding1)

        sol1 = solver(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg,
                     lam_x0=guess2, lam_g0=guess3) #, p=p_nlp0)

    print '\n****************************************************************'
    print '               Solved (Nominal MHE) Iteration No. ' + repr(i)
    print '****************************************************************\n'

    # Build the solve time array for nominal MHE solve
    elapsedTime = time.time() - startTime
    solve_time_mhe = append(solve_time_mhe, elapsedTime)
    
    # Now we use the predicted solve to resolve the actual MHE
    # Step 1: Get the predicted measurements
    # Simulate model one step using the estimated values
    print '\n Generating predicted measurements... \n'
    if i == 0:  # If first solve then use the initial values
        Ik_next_p = I_next(x0=X_ini,
                           p=vertcat(par_val_mismatch, U_dm2[i] + U_dm1[i] + w_dm[i]),
                           z0=Z_ini)
    else:
        Ik_next_p = I_next(x0=x_guess[:, -1],
                           p=vertcat(par_val_mismatch, U_dm2[i] + U_dm1[i] + w_dm[i]),
                           z0=z_guess[:, -1])

    # Get all the dependent variable
    d_all_next_p = f_y(Ik_next_p['xf'], Ik_next_p['zf'], U_dm2[i] + U_dm1[i] + w_dm[i],
                       par_val_mismatch)

    # Choose to get the measurement vector
    y_next_p = vertcat(d_all_next_p[318:321], d_all_next_p[322],
                       Ik_next_p['xf'][26:28],                  # continuous measurements
                       Ik_next_p['xf'][19], d_all_next_p[396],  # molten metal measurements
                       d_all_next_p[196], d_all_next_p[201], d_all_next_p[203],
                       d_all_next_p[205:207])                   # slag measurements

    # Assign the value to the measurement parameters
    # p_nlp0 = []  # Parameter assignment
    # p_nlp0 += y_mat_wn[0:Nmeas_end, i + 1].nonzeros()
    p_nlp0 = y_next_p[0:Nmeas_end] # .nonzeros()

    # Print statements
    print '\n****************************************************************'
    print '           Solving predicted problem, Iteration No. ' + repr(i)
    print '****************************************************************\n'
    # Step 2: Solve the predicted problem
    # Start clock to measure simulation time
    startTime = time.time()

    if i <= mhe_horizon-1:
        if i == 0:  # if first solve no dual guess available
            sol = solver(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg) #, p=p_nlp0)
        else:
            sol = solver(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg,
                         lam_x0=guess2, lam_g0=guess3) #, p=p_nlp0)
    else:  # nonbatch mhe
        sol = solver(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg,
                     lam_x0=guess2, lam_g0=guess3) #, p=p_nlp0)

    # Build the solve time array for predicted MHE solve
    elapsedTime = time.time() - startTime

    print '\n****************************************************************'
    print '              Solved (Predicted MHE) Iteration No. ' + repr(i)
    print '****************************************************************\n'

    solve_time_pmhe = append(solve_time_pmhe, elapsedTime)

    # Extract the predicted state estimated and the corresponding algebraic 
    # variable values
				
    # Primal solution extraction
    w_opt_p = sol['x'].full().flatten()	

    #%%# Extraction process from predicted results
    # First get the predicted estimates 
    x_est_p = w_opt_p[-(Nstates+Nalgvars):-(Nalgvars)]

    # Now get the algebraic variables :(																						
    # loop over the finite elements
    if i >= mhe_horizon-1:# non-batch mode
    
        # Drop the first element values
        w_opt_p = w_opt_p[(Nstates)*(steps+1)+(Nalgvars)*steps+Ninputs:]       
        
        # Repeating part
        w_opt_p1 = w_opt_p[Nstates + Nalgvars:]
        zp_guess = []
        
        for iw in range(mhe_horizon-1):
            # Select the first set of values for a finite element
            w_opt_p2 = DM(w_opt_p1[(Ninputs+(Nstates)*(steps+1)+(Nalgvars)*steps)*iw:\
                            (Ninputs+(Nstates)*(steps+1)+(Nalgvars)*steps)*(iw+1)])            
            
            # Remove the input part
            w_opt_p3 = w_opt_p2[Ninputs+Nstates:]
            
            for ix in range(steps):
                xz_hold = w_opt_p3[(Nstates+Nalgvars)*ix:(Nstates+Nalgvars)*(ix+1)]
                # x_guess = horzcat(x_guess,xz_hold[0:Nstates])
                zp_guess = horzcat(z_guess,xz_hold[Nstates:Nstates+Nalgvars])
                
    else: # batch model 
        # Repeating part
        w_opt_p1 = w_opt_p[Nstates + Nalgvars:]
        zp_guess = []
                     
        for iw in range(i+1):
            # Select the first set of values for a finite element
            w_opt_p2 = DM(w_opt_p1[(Ninputs+(Nstates)*(steps+1)+(Nalgvars)*steps)*iw:\
                            (Ninputs+(Nstates)*(steps+1)+(Nalgvars)*steps)*(iw+1)])            
            
            # Remove the input part
            w_opt_p3 = w_opt_p2[Ninputs+Nstates:]
            
            for ix in range(steps):
                xz_hold = w_opt_p3[(Nstates+Nalgvars)*ix:(Nstates+Nalgvars)*(ix+1)]
                # x_guess = horzcat(x_guess,xz_hold[0:Nstates])
                zp_guess = horzcat(zp_guess,xz_hold[Nstates:Nstates+Nalgvars]) 

    # Now get the last value from z_guess!
    z_est_p = zp_guess[:,-1]
																
    # Step 3: Resolve the original problem			
    p_nlp0 = y_mat_wn[0:Nmeas_end, i + 1] #.nonzeros()

    print '\n****************************************************************'
    print '          Solving actual fast problem, Iteration No. ' + repr(i)
    print '****************************************************************\n'

    # Start clock to measure simulation time
    startTime = time.time()

    # Resolve the actual problem with warm start
    # Warm start options for ipopt
    if i == 0:  # For other i these options will get defined from above
        opts_nlp["ipopt.warm_start_init_point"] = 'yes'
        opts_nlp["ipopt.warm_start_bound_push"] = 1e-6
        opts_nlp["ipopt.warm_start_mult_bound_push"] = 1e-6
        opts_nlp["ipopt.mu_init"] = 1e-6
        solver = nlpsol('solver', 'ipopt', prob, opts_nlp)

    sol = solver(x0=sol['x'], lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg,
                 lam_x0=sol['lam_x'], lam_g0=sol['lam_g']) #, p=p_nlp0)

    # Build the solve time array for nominal MHE solve
    elapsedTime = time.time() - startTime

    print '\n****************************************************************'
    print '                Solved (Fast MHE) Iteration No. ' + repr(i)
    print '****************************************************************\n'

    solve_time_fmhe = append(solve_time_fmhe, elapsedTime)

    # Primal solution extraction
    w_opt = sol1['x'].full().flatten()

    # Dual solution for x and g
    # We will extract what is needed and add zeros to fill the rest
    w_opt_dual = sol1['lam_x'].full().flatten()
    w_opt_lam_g = sol1['lam_g'].full().flatten()

    # Build the ipopt iteration count array
    solve_iter_mhe = append(solve_iter_mhe,
                            solver.stats().get('iter_count'))
    
    ####
    # Resolve with mhe_obj
    #prob = {'f': mhe_obj, 'x': vertcat(*w), 'g': vertcat(*g)}
    #solver = nlpsol('solver', 'ipopt', prob, opts_nlp)
    #sol = solver(x0=w_opt, lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg)
    
    #%%#################################################### 
    #          Post prosessing after MHE solve
    #######################################################
    
    ## Step 1: Update the arrival cost covariance using EKF update (non-batch MHE)           
    if i >= mhe_horizon-1:# and flag==0:
    
        # linearize around initial state values
        #f_c = dae.create('f_c', ['x','z','u','p'], ['jac_ddef_x','jac_ddef_z',
        #         'jac_alg_x','jac_alg_z','jac_ode_x','jac_ode_z'])
        my_e, my_f, my_g, my_h, my_i, my_j = f_c(w_opt[0:Nstates], 
                                           w_opt[Nstates:Nstates + Nalgvars],
                                           w_opt[Nstates + Nalgvars:Nstates + Nalgvars + Nstates + Ninputs],
                                           #U_dm2[i] + U_dm1[i] + w_dm[i],
                                           par_val_mismatch)
        # """                                   
        # Forward Euler Discretization method
        my_i = my_i + eye(Nstates)
        my_k = solve(my_h,my_g)
        my_a = my_i - mtimes(my_j,my_k)        

        # No need for exponentials
        a_c2d = my_a
        
        # get the full C matrix
        my_c1 = my_e - mtimes(my_f,my_k) 

        # Get the no of measurements at that time instant
        Nmeas = meas_struct[i - mhe_horizon + 1]
        
        # Build appropriate C matrix
        if Nmeas == 13:                   
            my_c = vertcat(my_c1[318,:], my_c1[319,:], my_c1[320,:], my_c1[322,:], 
                           my_c_t_roof, my_c_t_wall,
                           my_c_mm_t, my_c1[396,:],
                           my_c1[196,:], my_c1[201,:], my_c1[203,:], my_c1[205,:],
                           my_c1[206,:])
            
            # Apply EKF update formula (no inverse!)
            P = Q1 + mtimes(mtimes(a_c2d,
                           (P - mtimes(mtimes(mtimes(P,my_c.T),
                                              solve(R1_full + mtimes(mtimes(my_c,P),my_c.T), my_c)),P))),
                           a_c2d.T)
            S = linalg.inv(P)                
                           
        if Nmeas == 8:
            my_c = vertcat(my_c1[318,:], my_c1[319,:], my_c1[320,:], my_c1[322,:], 
                           my_c_t_roof, my_c_t_wall,
                           my_c_mm_t, my_c1[396,:])
                           
            # Apply EKF update formula (no inverse!)
            P = Q1 + mtimes(mtimes(a_c2d,
                           (P - mtimes(mtimes(mtimes(P,my_c.T),
                                              solve(R1_half + mtimes(mtimes(my_c,P),my_c.T), my_c)),P))),
                           a_c2d.T)
            S = linalg.inv(P)               
                           
        if Nmeas == 6:
            my_c = vertcat(my_c1[318,:], my_c1[319,:], my_c1[320,:], my_c1[322,:], 
                           my_c_t_roof, my_c_t_wall)
                           
            # Apply EKF update formula (no inverse!)
            P = Q1 + mtimes(mtimes(a_c2d,
                           (P - mtimes(mtimes(mtimes(P,my_c.T),
                                              solve(R1 + mtimes(mtimes(my_c,P),my_c.T), my_c)),P))),
                           a_c2d.T) 
            S = linalg.inv(P)                                               
    
    #%%# Extraction process from results
    # First get the estimates 
    x_mat_est = horzcat(x_mat_est, w_opt[-(Nstates+Nalgvars):-(Nalgvars)])

    # Now loop over the finite elements
    if i >= mhe_horizon-1:# non-batch mode
    
        # Drop the first element values
        w_opt = w_opt[(Nstates)*(steps+1)+(Nalgvars)*steps+Ninputs:]
        
        # Extract the initial state values
        X_ini = DM(w_opt[0:Nstates])
        
        # Extract the initial algebraic variable values
        Z_ini = DM(w_opt[Nstates:Nstates + Nalgvars])
        
        # Repeating part
        w_opt1 = w_opt[Nstates + Nalgvars:]
        u_opt = []
        x_guess = []
        z_guess = []
        
        for iw in range(mhe_horizon-1):
            # Select the first set of values for a finite element
            w_opt2 = DM(w_opt1[(Ninputs+(Nstates)*(steps+1)+(Nalgvars)*steps)*iw:\
                            (Ninputs+(Nstates)*(steps+1)+(Nalgvars)*steps)*(iw+1)])
            
            # Now extract the input values (model noise included :| )                
            u_opt = horzcat(u_opt,w_opt2[0:Ninputs+Nstates])
            
            # Remove the input part
            w_opt3 = w_opt2[Ninputs+Nstates:]
            
            for ix in range(steps):
                xz_hold = w_opt3[(Nstates+Nalgvars)*ix:(Nstates+Nalgvars)*(ix+1)]
                x_guess = horzcat(x_guess,xz_hold[0:Nstates])
                z_guess = horzcat(z_guess,xz_hold[Nstates:Nstates+Nalgvars])
                
    else: # batch model 
        X_ini = DM(w_opt[0:Nstates])
        
        # Extract the initial algebraic variable values
        Z_ini = DM(w_opt[Nstates:Nstates + Nalgvars])
        
        # Repeating part
        w_opt1 = w_opt[Nstates + Nalgvars:]
        x_guess = []
        z_guess = []
        u_opt = []
                     
        for iw in range(i+1):
            # Select the first set of values for a finite element
            w_opt2 = DM(w_opt1[(Ninputs+(Nstates)*(steps+1)+(Nalgvars)*steps)*iw:\
                            (Ninputs+(Nstates)*(steps+1)+(Nalgvars)*steps)*(iw+1)])
            
            # Now extract the input values (model noise included :| )               
            u_opt = horzcat(u_opt,w_opt2[0:Ninputs+Nstates])
            
            # Remove the input part
            w_opt3 = w_opt2[Ninputs+Nstates:]
            
            for ix in range(steps):
                xz_hold = w_opt3[(Nstates+Nalgvars)*ix:(Nstates+Nalgvars)*(ix+1)]
                x_guess = horzcat(x_guess,xz_hold[0:Nstates])
                z_guess = horzcat(z_guess,xz_hold[Nstates:Nstates+Nalgvars]) 
                
    # Check if solved successfully otherwise break
    if solver.stats().get('return_status')=='Solve_Succeeded':
        flag = 0
    else:
        flag = 1
        break

    #%%####################################################### 
    #                 Optimization Solve
    ##########################################################
    if i+1 in opt_call:  # If optimization called

        # Print statements
        print '\n****************************************************************'
        print '           Optimization called at Iteration No. '+ repr(i)
        print '****************************************************************\n'
        
        # Keep calculating the profit from injection to plant
    
        # Define mode:- 0: feasibility/ 1:optimization
        if i == 10:								
            oi_mode = 1
        else:												
            oi_mode  = 1
        
        # Specify the time where to start optimization
        oi_n_ini = i+1  
        
        # Specify relative time to start optimization
        oi_n_ini_adv = diff  # diff[count]
        count = count + 1        
        
        p_nlp1 = x_mat_est[:,-1]
        #p_nlp1 = x_mat_act[:,-1]							   
        oi_X_ini = p_nlp1  # x_mat_est[:,-1]
        oi_X_ini[-1] = 0  # Disturbance State
        oi_X_ini_be = p_nlp1  # x_mat_est[:,-1]
        oi_X_ini_be[-1] = 0  # Disturbance State
        
        p_nlp2 = z_guess[:,-1]
        #p_nlp2 = z_mat_act[:,-1]
        oi_Z_ini = p_nlp2  # z_guess[:,-1]                                          
        
        # No of control intervals
        oi_N = 60 - oi_n_ini
    
        #%%##########################################
        # Setup the optimization here
        #############################################
        
        # No of steps per control interval
        #if i in [10,23,24]:								
        if i in [7,  10,23,24,26]:
            oi_steps = 5 #(4 works fine),(1,2,3,5,6locally infeasible),(7:max iterations):'limited-memory'ipopt:N=60 
        # 21/4 :( 25/25/5 :( (after that only 4 & 8)
        else:
            if i in [21]: 
                oi_steps = 6
            else:																
                oi_steps = 4                
        # step size
        oi_h = 1/oi_steps
        
        # output grid
        time_grid = linspace(0, 1, oi_steps+1)
        
        # setup the idas integrator
        oi_I_full = integrator('integrator', 'idas', dae1,
                               {'tf':1, 'calc_ic':True,
							 'abstol':1e-6, 'reltol':1e-6, # earlier was 8 																																
                                'grid': time_grid, 'output_t0':True})
        
        # Start with an empty NLP
        w=[]
        w0 = []
        lbw = []
        ubw = []
        J = 0
        g=[]
        lbg = []
        ubg = []
        
        # "Lift" initial conditions
        oi_X_ini_def = SX.sym('oi_X_ini_def', Nstates)
        w += [oi_X_ini_def]
        lbw += [k for k in oi_X_ini.nonzeros()]
        ubw += [k for k in oi_X_ini.nonzeros()] 
        w0 += [k for k in oi_X_ini.nonzeros()]
        
        oi_Xk = oi_X_ini_def
        
        # Define contrainer for exact state values obtained from IDAS
        oi_X_exact = oi_X_ini
        oi_Z_exact = oi_Z_ini
    
        #%%############################################################
        # Back euler loop starts here
        ############################################################### 
        for oi_j1 in range(oi_N):
                
            # Initialization of backeuler simulation
												
            if oi_j1 != 70:												
		       # Apply newton solve to get consistent initial conditions
                x_1_new = s11(oi_Z_ini, oi_X_ini, 
                          vertcat(par_val_mismatch, opt_u_nom[:,oi_j1+oi_n_ini]))

                # Assign the new consistent alg var values 
                oi_Z_ini = x_1_new[0]
									
            oi_Ik1 = oi_I_full(x0=oi_X_ini, 
                         p=vertcat(par_val_mismatch,
                                   opt_u_nom[:,oi_j1+oi_n_ini]),
                                   z0 = oi_Z_ini)
                                   
            oi_X_ini = oi_Ik1['xf'][:,-1]
            oi_Z_ini = oi_Ik1['zf'][:,-1]
            oi_X_exact = horzcat(oi_X_exact, oi_X_ini)
            oi_Z_exact = horzcat(oi_Z_exact, oi_Z_ini)
            
            # NLP part
            # New NLP variable for the control
            oi_Uk = SX.sym('oi_U_' + str(oi_j1), Nstates + Ninputs)
            w   += [oi_Uk]
            
            # Fix inputs according to the feasibility/optimization
            if oi_mode == 0:
                # Inputs fixed as we want to solve for feasibility only
                lbw += U_dm2[oi_j1+oi_n_ini] + U_dm1[oi_j1+oi_n_ini]
                ubw += U_dm2[oi_j1+oi_n_ini] + U_dm1[oi_j1+oi_n_ini]
            else:    
                lbw += [x1_ * 1 for x1_ in U_dm2[oi_j1+oi_n_ini]]
                lbw += [0.8*U_dm1[oi_j1+oi_n_ini][0]]
                lbw += [0.7*U_dm1[oi_j1+oi_n_ini][1]]#ss_p_arc 0.8
                lbw += [x2_ * 0.8 for x2_ in U_dm1[oi_j1+oi_n_ini][2:]]
            
                ubw += [x3_ * 1 for x3_ in U_dm2[oi_j1+oi_n_ini]]
                ubw += [1.2*U_dm1[oi_j1+oi_n_ini][0]]  
                ubw += [1.3*U_dm1[oi_j1+oi_n_ini][1]] #ss_p_arc 1.2
                ubw += [x4_ * 1.2 for x4_ in U_dm1[oi_j1+oi_n_ini][2:]]
                
                # Add contribution of inputs to cost function
                J = J  + u7_coef*oi_Uk[7] + (u8_coef[oi_j1+oi_n_ini]*oi_Uk[8] + u9_coef*(oi_Uk[9]+oi_Uk[10]+oi_Uk[11]))
                #J = J  + 1.5*0.213650*Uk[7]/60 + (0.0585*1000*Uk[8]/60 + 1.5*0.07875*(Uk[9]+Uk[10]+Uk[11])/60)
    
            # Provide bounds for model noises 
            lbw += [0]*Nstates
            ubw += [0]*Nstates

            if i in [7,8, 10,11, 21,22, 23,25, 26,27]:
                w0  += U_dm2[o_j1+o_n_ini] + U_dm1[o_j1+o_n_ini] + w_dm[o_j1+o_n_ini]
            
            for oi_j in range(oi_steps):
        
                oi_X_ini_be = oi_Ik1['xf'][:,oi_j+1]
                oi_Z_ini_be = oi_Ik1['zf'][:,oi_j+1]
                
                # NLP part
                # State at the euler points
                oi_Xk_end = SX.sym('oi_X_' + str(oi_j1+1) + '_' + str(oi_j), Nstates)
                oi_Zk_end = SX.sym('oi_Zz_' + str(oi_j1+1) + '_' + str(oi_j), Nalgvars)
                
                w   += [oi_Xk_end]
        
                for i5 in range(Nstates):
                    if i5 in range(Nstates-1):#[26,27,7]:#,7,8,10]:
                        if abs(oi_X_ini_be[i5]) > 0 and abs(oi_X_ini_be[i5]) < 5:
                            if oi_X_ini_be[i5] < 0:
                                lbw += [-10]
                            else:
                                lbw += [0]
                        else:
                            if oi_X_ini_be[i5] <= 0:
                                lbw += (1.5*oi_X_ini_be[i5]).nonzeros()
                            else:    
                                lbw += (0.5*oi_X_ini_be[i5]).nonzeros()
                    else:  # 0 for disturbance state   
                        # lbw+= [X_ini_be.nonzeros()[i5]] 
                        lbw += [0]
                   
                for i4 in range(Nstates):
                    if i4 in range(Nstates-1):#[26,27,7]:#,7,8,10]:            
                        if abs(oi_X_ini_be[i4]) > 0 and abs(oi_X_ini_be[i4]) < 5:
                            if oi_X_ini_be[i4] < 0:
                                ubw += [0]
                            else:
                                ubw += [10]
                        else:
                            if oi_X_ini_be[i4] < 0:
                                ubw += (0.5*oi_X_ini_be[i4]).nonzeros()
                            else:    
                                ubw += (1.5*oi_X_ini_be[i4]+0.0001).nonzeros()  
                    else:  # 0 for disturbance state
                        # ubw += [X_ini_be.nonzeros()[i4]] 
                        ubw += [0]
                if i in [7,8, 10,11,21,22,23,25,26,27]:       
                    w0 += [k for k in oi_X_ini_be.nonzeros()]
                
                w   += [oi_Zk_end]
            
                for i3 in range(Nalgvars):
                    if i3 in range(14)+[15,17,18,20,22,24,26]+range(28,43)+range(52,108)+[110,113,115,119,120]: 
                        if abs(oi_Z_ini_be[i3]) > 0 and abs(oi_Z_ini_be[i3]) < 5:
                            if oi_Z_ini_be[i3] < 0:
                                lbw += [-10]
                            else:
                                lbw += [0]
                        else:
                            if oi_Z_ini_be[i3] <= 0:
                                lbw += (1.5*oi_Z_ini_be[i3]).nonzeros()
                            else:    
                                lbw += (0.5*oi_Z_ini_be[i3]).nonzeros()
                    #else:      #  0     1   2   3    4    5   6   7    8   9   10  11   12  13   14   15     
                    #    if i3 in [14*, 51, 16, 44*, 43*, 19, 45, 21*, 46, 23*, 47, 25*, 48, 27*, 49*, 50*]
                    else:      #  0     1   2   3    4    5   6   7    8   9   10  11   12  13   14   15     
                        if i3 in [51, 16, 19, 45, 46, 47, 48]:
                            lbw += (1.2*oi_Z_ini_be[i3]).nonzeros() # earlier 1.2
                        else:
                            lbw += (1.5*oi_Z_ini_be[i3]).nonzeros() # earlier 1.5
                            #lbw += [-1e10]
        
                for i2 in range(Nalgvars):
                    if i2 in range(14)+[15,17,18,20,22,24,26]+range(28,43)+range(52,108)+[110,113,115,119,120]: 
                        if abs(oi_Z_ini_be[i2]) > 0 and abs(oi_Z_ini_be[i2]) < 5:
                            if oi_Z_ini_be[i2] < 0:
                                ubw += [0]
                            else:
                                ubw += [10]
                        else:
                            if oi_Z_ini_be[i2] < 0:
                                ubw += (0.5*oi_Z_ini_be[i2]).nonzeros()
                            else:    
                                ubw += (1.5*oi_Z_ini_be[i2]+0.0001).nonzeros()
                    else:      #  0     1   2   3    4    5   6   7    8   9   10  11   12  13   14   15     
                        if i2 in [51, 16, 19, 45, 46, 47, 48]:
                            ubw += (0.8*oi_Z_ini_be[i2]).nonzeros()#earlier 0.8
                        else:
                            #ubw += (0.9*Z_ini_be[i2]).nonzeros()                    
                            #ubw += [0.001] # [0.0001] upto 47      
                            ubw += [0.0001]# upto 47                       
                if i in [7,8, 10,11,21,22,23,25,26,27]:                            
                    w0  += [k for k in oi_Z_ini_be.nonzeros()]
            
                # Euler equations
                #f = Function('f', [x_all, z_all, params], [ode, alg])
                oi_fj, oi_zj = f(oi_Xk_end, oi_Zk_end, par_val_mismatch, oi_Uk)
                    
                g += [oi_Xk[0:(Nstates-1)] - oi_Xk_end[0:(Nstates-1)] + oi_fj[0:(Nstates-1)]*oi_h]
                
                lbg += [0]*(Nstates-1)
                ubg += [0]*(Nstates-1)
                
                g += [oi_zj]
                lbg += [0]*Nalgvars
                ubg += [0]*Nalgvars            
            
                # Add equality constraint
                oi_Xk = oi_Xk_end
        
        if oi_mode == 1: # if mode is optimization then add objective contribution
            #J = J - 0.359472*Xk_end[12]
    
            J = J - xk15o_coef*((oi_Xk_end[15]*xk15i_coef*xk15g_coef) + (((oi_Zk_end[42])*exp(oi_Zk_end[44]))*xk15i_coef))
    
            #J = J - 0.359472*((o_Xk_end[15]*55.85*0.1) + (((o_Zk_end[42])*exp(o_Zk_end[44]))*55.85))
        #%%#################################################
        # Solve the NLP
        ###################################################    
        
        # NLP solver options
        opts_nlp = {}
        opts_nlp["ipopt.linear_solver"] = 'ma27'
        #opts_nlp["ipopt.linear_system_scaling"] = 'slack-based'
        #opts_nlp["ipopt.alpha_for_y"]='primal-and-full'
        #opts_nlp["ipopt.obj_scaling_factor"] = 1e-3
        opts_nlp["ipopt.hessian_approximation"] = 'limited-memory'
        opts_nlp["ipopt.max_iter"] = 300								
                
        opts_nlp["ipopt.warm_start_init_point"] = 'yes'
        opts_nlp["ipopt.warm_start_bound_push"] = 1e-6
        opts_nlp["ipopt.warm_start_mult_bound_push"] = 1e-6
        #opts_nlp["ipopt.warm_start_slack_bound_push"] = 1e-6
        opts_nlp["ipopt.mu_init"] = 1e-6 
        #opts_nlp["ipopt.warm_start_mult_init_max"] = 10
            
        prob = {'f': J, 'x': vertcat(*w), 'g': vertcat(*g)}
        solver = nlpsol('solver', 'ipopt', prob, opts_nlp);        
        # i =21,23,24 problem, check 26
        #if i in range(10) + range(12,21) + [24,26,27,28,29]:        
        #if i in range(10) + range(12,21) + [24]+ range(26,35):  
        if i in range(7) + [9] + range(12,21) + [24]+ range(28,n_mhe-1): 
            w0 = append(w0, 
                    o_w_opt_nom1[(Nstates+Nalgvars)*(oi_steps*oi_n_ini_adv)\
                    +(Nstates+Ninputs)*oi_n_ini_adv:])
        
            lamx0 = []
            lamg0 = []
            lamx0 += [0]*Nstates 
								
            lamx0 = append(lamx0, opt_sol_nom['lam_x'][(Nstates+Nalgvars)*(oi_steps*oi_n_ini_adv)\
                        +(Nstates+Ninputs)*oi_n_ini_adv + Nstates:])
            lamg0 = opt_sol_nom['lam_g'][(Nstates+Nalgvars-1)*(oi_steps*oi_n_ini_adv):]            
        
            # Start clock to measure simulation time
            startTime = time.time()
        
            # Solve the NLP in warmstart mode
            if i in range(6,n_mhe-1):  # range(6,11)
            #if i in range(6,30):  # range(6,11)													
                opt_sol_1 = solver(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, 
                             ubg=ubg, lam_x0=lamx0, lam_g0=lamg0) 
            else:										
                opt_sol_1 = solver(x0=w0, lbx=append(x_mat_est[:,-1], o_lbw[(Nstates+Nalgvars)*(oi_steps*count)\
                        +(Nstates+Ninputs)*count + Nstates:]),
						    ubx=append(x_mat_est[:,-1], o_ubw[(Nstates+Nalgvars)*(oi_steps*count)\
                        +(Nstates+Ninputs)*count + Nstates:]),
                               lbg=lbg, ubg=ubg, lam_x0=lamx0, lam_g0=lamg0)                     
        
            # Stop the clock and read in the total time taken                    
            elapsedTime = time.time() - startTime
				
        else:
            # If MPC called with changed finite elements or at the next time 
            # step then no previous solve info used									
            startTime = time.time()									
            opt_sol_1 = solver(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, 
                             ubg=ubg)	
            elapsedTime = time.time() - startTime

        # Build the solve time array
        solve_time_opt = append(solve_time_opt, elapsedTime)
												
        print '\n****************************************************************'
        print '               Solved (Nominal MPC) Iteration No. ' + repr(i)
        print '****************************************************************\n'												

        print '\n Geting predicted state estimates... \n' 								
        p_nlp1 = x_est_p
        p_nlp2 = z_est_p 
								
        # Print statements
        print '\n****************************************************************'
        print '       Solving predicted MPC problem, Iteration No. ' + repr(i)
        print '****************************************************************\n'
        # Step 2: Solve the predicted problem
        # Start clock to measure simulation time
        #if i in range(10) + range(12,21) + [24]+ range(28,n_mhe-1):
        if i in range(7) + [9] + range(12,21) + [24]+ range(28,n_mhe-1):    
        #if i in range(10) + range(12,21) + [24,26]: 
            startTime = time.time()

            # Solve the NLP in warmstart mode
            if i in range(6,n_mhe-1):				
                opt_sol_2 = solver(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, 
                             ubg=ubg, lam_x0=lamx0, lam_g0=lamg0) 
            else:										
                opt_sol_2 = solver(x0=w0, lbx=append(p_nlp1, o_lbw[(Nstates+Nalgvars)*(oi_steps*count)\
                        +(Nstates+Ninputs)*count + Nstates:]),
						    ubx=append(p_nlp1, o_ubw[(Nstates+Nalgvars)*(oi_steps*count)\
                        +(Nstates+Ninputs)*count + Nstates:]),
                               lbg=lbg, ubg=ubg, lam_x0=lamx0, lam_g0=lamg0)  

            # Build the solve time array for predicted MHE solve
            elapsedTime = time.time() - startTime
        else:
            startTime = time.time()		   						
            opt_sol_2 = solver(x0=w0, lbx=lbw, ubx=ubw, lbg=lbg, 
                             ubg=ubg)
            elapsedTime = time.time() - startTime
												
        print '\n****************************************************************'
        print '              Solved (Predicted MPC) Iteration No. ' + repr(i)
        print '****************************************************************\n'

        solve_time_pmpc = append(solve_time_pmpc, elapsedTime)

        # Step 3: Resolve the original problem
        #p_nlp1 = x_mat_est[:,-1]
        #p_nlp2 = z_guess[:,-1] 

        print '\n****************************************************************'
        print '           Solving actual fast MPC, Iteration No. ' + repr(i)
        print '****************************************************************\n'

        # Start clock to measure simulation time
        startTime = time.time()

        # Resolve the actual problem with warm start
        #if i in range(10) + range(12,21) + [24]+ range(28,n_mhe-1):
        if i in range(7) + [9] + range(12,21) + [24]+ range(28,n_mhe-1):
            if i in range(6,n_mhe-1):				
                opt_sol_nom = solver(x0=opt_sol_2['x'], lbx=lbw, ubx=ubw, lbg=lbg, 
                              ubg=ubg, lam_x0=opt_sol_2['lam_x'], 
                              lam_g0=opt_sol_2['lam_g']) 
            else:										
                opt_sol_nom = solver(x0=opt_sol_2['x'], lbx=append(p_nlp1, o_lbw[(Nstates+Nalgvars)*(oi_steps*count)\
                        +(Nstates+Ninputs)*count + Nstates:]),
						    ubx=append(p_nlp1, o_ubw[(Nstates+Nalgvars)*(oi_steps*count)\
                        +(Nstates+Ninputs)*count + Nstates:]),
                               lbg=lbg, ubg=ubg, lam_x0=opt_sol_2['lam_x'], 
                               lam_g0=opt_sol_2['lam_g'])
        else:	
            opt_sol_nom = solver(x0=opt_sol_2['x'], lbx=lbw, ubx=ubw, lbg=lbg, 
                              ubg=ubg, lam_x0=opt_sol_2['lam_x'], 
                              lam_g0=opt_sol_2['lam_g']) 											
        # Build the solve time array for nominal MPC solve
        elapsedTime = time.time() - startTime

        print '\n****************************************************************'
        print '                Solved (Fast MPC) Iteration No. ' + repr(i)
        print '****************************************************************\n'

        solve_time_fmpc = append(solve_time_fmpc, elapsedTime)								
								
        opt_sol_nom_f = opt_sol_nom['x'].full().flatten() 
        o_w_opt_nom1 = opt_sol_nom_f[Nstates:]               

        # Extract the optimal inputs from the solution and replace in opt_u_nom
        # Repeating part        
        opt_u_nom_temp = []
                             
        for iw in range(oi_N):
            # Select the first set of values for a finite element
            o_w_opt_nom2 = DM(o_w_opt_nom1[(Ninputs+(Nstates)*(oi_steps+1)+(Nalgvars)*oi_steps)*iw:\
                               (Ninputs+(Nstates)*(oi_steps+1)+(Nalgvars)*oi_steps)*(iw+1)])
                    
            # Now extract the input values (model noise included :| )               
            opt_u_nom_temp = horzcat(opt_u_nom_temp, o_w_opt_nom2[0:Ninputs+Nstates])           
 
        # Replace the temp in opt_u_nom
        opt_u_nom[:,oi_n_ini:] = opt_u_nom_temp

###################################################
#              Cost calculations
###################################################
# Add last input cost to plant objective function
J_cl = J_cl  + u7_coef*opt_u_nom[7,59] + \
          (u8_coef2[59]*opt_u_nom[8,59] + u9_coef*(opt_u_nom[9,59]+opt_u_nom[10,59]+opt_u_nom[11,59]))

J_cl_arc_power = J_cl_arc_power + u8_coef2[59]*opt_u_nom[8,59]
J_cl_other = J_cl - J_cl_arc_power

# Add profit to plant objective function
# Simulate the plant one step 
Ik_next = I_next(x0 = x_ini_start_mhe, 
                 p = vertcat(par_val, opt_u_nom[:,59]),
                 z0 = z_ini_start_mhe)  
 
J_cl_earnings = xk15o_coef*((Ik_next['xf'][15]*xk15i_coef*xk15g_coef) + (((Ik_next['zf'][42])*exp(Ik_next['zf'][44]))*xk15i_coef))                              
J_cl = J_cl - xk15o_coef*((Ik_next['xf'][15]*xk15i_coef*xk15g_coef) + (((Ik_next['zf'][42])*exp(Ik_next['zf'][44]))*xk15i_coef))

# Add the last measurement to y_mat
# Get all the dependent variable values
d_all_end = f_y(Ik_next['xf'][:,0], Ik_next['zf'][:,0], opt_u_nom[:,59], 
            par_val)
                  
# Choose to get the measurement vector and add noise accordingly
y_end = vertcat(d_all_end[318:321], d_all_end[322], Ik_next['xf'][26:28,0],# continuous measurements
                  Ik_next['xf'][19,0], d_all_end[396], # molten metal measurements 
                  d_all_end[196], d_all_end[201], d_all_end[203], d_all_end[205:207]) #slag measurements          

																    # Build the measurement matrix
y_mat = horzcat(y_mat, y_end)
###################################################
#                Plotting
################################################### 

# Define time grid
tgrid = range(x_mat_est.shape[1]) 

matplotlib.rcParams.update({'font.size': 8})
fig = plt.figure(1)
fig.set_size_inches(10.5, 12.5)

## Decide the subplot  
sp_c = 4
if disturb_add ==1:
    sp_r = (x_mat_est.shape[0]-1)/sp_c + 1
else:    
    sp_r = (x_mat_est.shape[0]-0)/sp_c# + 1
    
for i in range(Nstates):
    fig.add_subplot(sp_r,sp_c,i+1).plot(tgrid, x_mat_est[i,:].full().flatten(), 'b-')#,linewidth=4.0)
    fig.add_subplot(sp_r,sp_c,i+1).plot(tgrid, x_mat_act[i,:x_mat_est.shape[1]].full().flatten(), 'r--')
    xlabel('t')
    ylabel(dae.x[i].name())
    #legend(['est','tru'])

plt.tight_layout()
fig = plt.gcf()
fig.savefig('estimate_jpc.png')

# PLot the solve time
bar_count = len(solve_time_mhe)
fig1 = plt.figure(2)
fig1.set_size_inches(9, 2.5)
axis = range(bar_count)
width = 1/2.5
plt.bar(axis, solve_time_fmhe, width, color="blue", align='center')
xlabel('Solve iteration')
ylabel('CPU time (sec)')

fig1 = plt.gcf()
fig1.savefig('solvetime_mhe_jpc.png')

# PLot the MPC solve time
bar_count1 = len(solve_time_opt)
fig1 = plt.figure(3)
fig1.set_size_inches(9, 2.5)
axis = range(bar_count1)
width = 1/2.5
plt.bar(axis, solve_time_opt, width, color="blue", align='center')
xlabel('MPC Solve iteration')
ylabel('CPU time (sec)')

fig2 = plt.gcf()
fig2.savefig('solvetime_mpc_jpc.png')
# PLot the solve iteration count
#fig2 = plt.figure(3)
#fig2.set_size_inches(9, 2.5)
#plt.bar(axis, solve_iter_mhe, width, color="red", align='center')
#xlabel('Solve iteration')
#ylabel('No. of iterations')

#fig2 = plt.gcf()
#fig2.savefig('solveiter_new3_base.png')

plt.show()               
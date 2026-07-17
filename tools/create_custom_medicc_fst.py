# %%
import medicc
import fstlib
from medicc.factory import _get_int_cns_from_symbol_table


separator = "X"
max_cn = 8
max_pre_wgd_losses = 7

# FST weights
w_stay_gain = 0
w_open_gain = 40
w_extend_gain = 4

w_stay_loss = 0
w_open_loss = w_open_gain
w_extend_loss = w_extend_gain

wgd_cost = w_open_gain

SYMBOL_TABLE = medicc.create_symbol_table(max_cn=max_cn, separator=separator)
n = len(medicc.factory._get_int_cns_from_symbol_table(SYMBOL_TABLE, separator))

# %% 
G1step_chrom = medicc.factory.create_1step_gain_fst_whole_chrom(SYMBOL_TABLE, separator, w_open_gain, minimize=False)
G_chrom = medicc.create_nstep_fst(n - 1, G1step_chrom)

# %% 

L_LOH_1step_chrom = medicc.factory.create_1step_loss_fst_whole_chrom(SYMBOL_TABLE, separator, w_open_loss, minimize=False)
L_LOH_chrom = medicc.create_nstep_fst(n - 1, L_LOH_1step_chrom)

# %%
G1step = ~medicc.factory.create_1step_del_fst(SYMBOL_TABLE,
                                              separator=separator,
                                              exclude_zero=True,
                                              w_stay=w_stay_gain,
                                              w_open=w_open_gain,
                                              w_extend=w_extend_gain)

L_LOH_1step = medicc.factory.create_1step_del_fst(SYMBOL_TABLE,
                                                  separator=separator,
                                                  exclude_zero=False,
                                                  w_stay=w_stay_loss,
                                                  w_open=w_open_loss,
                                                  w_extend=w_extend_loss)

G = medicc.create_nstep_fst(n - 1, G1step)
L_LOH = medicc.create_nstep_fst(n - 1, L_LOH_1step)

# %%
# XX = fstlib.encode_determinize_minimize(L_LOH * G) # NOT gonna work
XX = L_LOH * G


# %%
W_1step = medicc.factory.create_1step_WGD_fst(SYMBOL_TABLE, separator,
                                              wgd_cost=wgd_cost,
                                              minimize=False,
                                              wgd_x2=False,
                                              total_cn=False)

MEDICC2_FST_with_extend_modeling = W_1step * G_chrom * L_LOH_chrom * XX


# %%
MEDICC2_FST_with_extend_modeling.write(f"/projects/ag-schwarzr/project-medicc/medicc2/medicc/objects/gain_loss_extend_open_{w_open_gain}_extend_{w_extend_gain}.fst")
XX.write(f"/projects/ag-schwarzr/project-medicc/medicc2/medicc/objects/gain_loss_extend_open_{w_open_gain}_extend_{w_extend_gain}_no_WGD.fst")
# %%

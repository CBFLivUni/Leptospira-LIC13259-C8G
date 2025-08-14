hide lines, all 
show cartoon

bg_color white
set cartoon_fancy_helices, 1

select chain A; color grey10, sele
select chain B; color grey90, sele

select chainB_rank0, chain B and ranked_0_Leptospira_borgpetersenii_LIC13259_no_his_and_Bos_taurus_C8_A8YXZ2
select chainB_rank1, chain B and ranked_1_Leptospira_borgpetersenii_LIC13259_no_his_and_Bos_taurus_C8_A8YXZ2
select chainB_rank2, chain B and ranked_2_Leptospira_borgpetersenii_LIC13259_no_his_and_Bos_taurus_C8_A8YXZ2

align chainB_rank1, chainB_rank0
align chainB_rank2, chainB_rank0

##############

select chainB_human, chain B and Leptospira_interrogans_s.Copenhageni_LIC13259_no_his_and_Homo_sapiens_C8_P07360
select chainB_bovine, chain B and Leptospira_interrogans_s.Copenhageni_LIC13259_no_his_and_Bos_taurus_C8_A8YXZ2

align chainB_human, chainB_bovine

#################

select cluster19, chain B and ranked_21_Leptospira_interrogans_LIC13259_no_his_and_Bos_taurus_C8_A8YXZ2
select cluster35, chain B and ranked_3_Leptospira_borgpetersenii_LIC13259_no_his_and_Bos_taurus_C8_A8YXZ2
select cluster36, chain B and ranked_3_Leptospira_interrogans_LIC13259_no_his_and_Bos_taurus_C8_A8YXZ2
select Cluster38, chain B and ranked_4_Leptospira_biflexa_LIC13259_no_his_and_Bos_taurus_C8_A8YXZ2
select cluster43, chain B and ranked_7_Leptospira_borgpetersenii_LIC13259_no_his_and_Bos_taurus_C8_A8YXZ2
select cluster49, chain B and ranked_9_Leptospira_interrogans_LIC13259_no_his_and_Bos_taurus_C8_A8YXZ2

align cluster19, cluster43
align cluster35, cluster43
align cluster36, cluster43
align cluster38, cluster43
align cluster49, cluster43


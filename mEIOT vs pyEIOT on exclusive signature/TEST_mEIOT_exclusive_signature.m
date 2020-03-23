function TEST_mEIOT_exclusive_signature

load('RESULTS_B100_EIOT_vs_others_20200312.mat')

[eiot_obj]=eiot_build(snv(SPEC_CAL_OLD),Y_CAL_OLD_1)
eiot_obj.S_I=RES_mean;
eiot_obj.S_E=[eiot_obj.S_hat;eiot_obj.S_I]
eiot_obj.num_si=11;
eiot_obj.num_e_si=11;
for i=1:178
[r_hat_sal(i,:),ri_hat_sal(i,:),ssr,m] = eiot_calc(snv(SPEC_B100_COMBINE(i,:)),eiot_obj);
end
imagesc(1:11,1:178,ri_hat_sal)
ylabel('sample #')
xlabel('batch #')
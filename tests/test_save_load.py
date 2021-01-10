import numpy as np
import fathon
from fathon import fathonUtils as fu
import os

# FUNCTIONALITY TESTS
# -------------------

TESTS_PATH = os.path.dirname(os.path.abspath(__file__))

#####
# Time series to test functionality
#####
ts = np.random.randn(1000)
#####

def get_object_path(name, ext=False):
    if ext:
        return os.path.join(TESTS_PATH, name + '.fathon')
    else:
        return os.path.join(TESTS_PATH, name)

#####
# Functionality test 1
# Save DFA object and reload
#####
def test_dfa_save_load():
    dfa = fathon.DFA(fu.toAggregated(ts))
    # save and load with empty results
    dfa.saveObject(get_object_path('dfa_obj'))
    n_load = fu.getObjectMember(get_object_path('dfa_obj', ext=True), 'n')
    F_load = fu.getObjectMember(get_object_path('dfa_obj', ext=True), 'F')
    assert np.array_equal(n_load, [])
    assert np.array_equal(F_load, [])
    #save and load with results
    n, F = dfa.computeFlucVec(fu.linRangeByStep(10, 500))
    H, I = dfa.fitFlucVec(100, 300)
    dfa.saveObject(get_object_path('dfa_obj'))
    n_load = fu.getObjectMember(get_object_path('dfa_obj', ext=True), 'n')
    F_load = fu.getObjectMember(get_object_path('dfa_obj', ext=True), 'F')
    assert np.array_equal(n_load, n)
    assert np.array_equal(F_load, F)
    # DFA from file
    dfa_2 = fathon.DFA(get_object_path('dfa_obj', ext=True))
    H_2, I_2 = dfa_2.fitFlucVec(100, 300)
    assert H_2 == H
    assert I_2 == I

#####
# Functionality test 2
# Save MFDFA object and reload
#####
def test_mfdfa_save_load():
    mfdfa = fathon.MFDFA(fu.toAggregated(ts))
    # save and load with empty results
    mfdfa.saveObject(get_object_path('mfdfa_obj'))
    n_load = fu.getObjectMember(get_object_path('mfdfa_obj', ext=True), 'n')
    F_load = fu.getObjectMember(get_object_path('mfdfa_obj', ext=True), 'F')
    q_list_load = fu.getObjectMember(get_object_path('mfdfa_obj', ext=True), 'qList')
    list_h_load = fu.getObjectMember(get_object_path('mfdfa_obj', ext=True), 'listH')
    assert np.array_equal(n_load, [])
    assert np.array_equal(F_load, [])
    assert np.array_equal(q_list_load, [])
    assert np.array_equal(list_h_load, [])
    #save and load with results
    n, F = mfdfa.computeFlucVec(fu.linRangeByStep(10, 500), fu.linRangeByStep(-1, 1))
    H, I = mfdfa.fitFlucVec(100, 300)
    mfdfa.saveObject(get_object_path('mfdfa_obj'))
    n_load = fu.getObjectMember(get_object_path('mfdfa_obj', ext=True), 'n')
    F_load = fu.getObjectMember(get_object_path('mfdfa_obj', ext=True), 'F')
    q_list_load = fu.getObjectMember(get_object_path('mfdfa_obj', ext=True), 'qList')
    list_h_load = fu.getObjectMember(get_object_path('mfdfa_obj', ext=True), 'listH')
    assert np.array_equal(n_load, n)
    assert np.array_equal(F_load, F)
    assert np.array_equal(q_list_load, fu.linRangeByStep(-1, 1))
    assert np.array_equal(list_h_load, H)
    # MFDFA from file
    mfdfa_2 = fathon.MFDFA(get_object_path('mfdfa_obj', ext=True))
    H_2, I_2 = mfdfa_2.fitFlucVec(100, 300)
    assert np.array_equal(H_2, H)
    assert np.array_equal(I_2, I)

#####
# Functionality test 3
# Save HT object and reload
#####
def test_ht_save_load():
    ht = fathon.HT(fu.toAggregated(ts))
    # save and load with empty results
    ht.saveObject(get_object_path('ht_obj'))
    ht_load = fu.getObjectMember(get_object_path('ht_obj', ext=True), 'ht')
    assert np.array_equal(ht_load, [])
    #save and load with results
    ht_mtx = ht.computeHt(np.array([100, 300], dtype=np.int64))
    ht.saveObject(get_object_path('ht_obj'))
    ht_load = fu.getObjectMember(get_object_path('ht_obj', ext=True), 'ht')
    assert np.array_equal(ht_load, ht_mtx)

#####
# Functionality test 4
# Save DCCA object and reload
#####
def test_dcca_save_load():
    dcca = fathon.DCCA(fu.toAggregated(ts), fu.toAggregated(ts))
    # save and load with empty results
    dcca.saveObject(get_object_path('dcca_obj'))
    n_load = fu.getObjectMember(get_object_path('dcca_obj', ext=True), 'n')
    F_load = fu.getObjectMember(get_object_path('dcca_obj', ext=True), 'F')
    n_rho_load = fu.getObjectMember(get_object_path('dcca_obj', ext=True), 'nRho')
    rho_load = fu.getObjectMember(get_object_path('dcca_obj', ext=True), 'rho')
    n_thr_load = fu.getObjectMember(get_object_path('dcca_obj', ext=True), 'nThr')
    conf_up_load = fu.getObjectMember(get_object_path('dcca_obj', ext=True), 'confUp')
    conf_down_load = fu.getObjectMember(get_object_path('dcca_obj', ext=True), 'confDown')
    assert np.array_equal(n_load, [])
    assert np.array_equal(F_load, [])
    assert np.array_equal(n_rho_load, [])
    assert np.array_equal(rho_load, [])
    assert np.array_equal(n_thr_load, [])
    assert np.array_equal(conf_up_load, [])
    assert np.array_equal(conf_down_load, [])
    #save and load with results
    n, F = dcca.computeFlucVec(fu.linRangeByStep(10, 200))
    H, I = dcca.fitFlucVec(100, 150)
    dcca.saveObject(get_object_path('dcca_obj'))
    n_load = fu.getObjectMember(get_object_path('dcca_obj', ext=True), 'n')
    F_load = fu.getObjectMember(get_object_path('dcca_obj', ext=True), 'F')
    n_rho_load = fu.getObjectMember(get_object_path('dcca_obj', ext=True), 'nRho')
    rho_load = fu.getObjectMember(get_object_path('dcca_obj', ext=True), 'rho')
    n_thr_load = fu.getObjectMember(get_object_path('dcca_obj', ext=True), 'nThr')
    conf_up_load = fu.getObjectMember(get_object_path('dcca_obj', ext=True), 'confUp')
    conf_down_load = fu.getObjectMember(get_object_path('dcca_obj', ext=True), 'confDown')
    assert np.array_equal(n_load, n)
    assert np.array_equal(F_load, F)
    assert np.array_equal(n_rho_load, [])
    assert np.array_equal(rho_load, [])
    assert np.array_equal(n_thr_load, [])
    assert np.array_equal(conf_up_load, [])
    assert np.array_equal(conf_down_load, [])
    # DCCA from file
    dcca_2 = fathon.DCCA(get_object_path('dcca_obj', ext=True))
    H_2, I_2 = dcca_2.fitFlucVec(100, 150)
    assert np.array_equal(H_2, H)
    assert np.array_equal(I_2, I)
    # rho
    n_rho, rho = dcca.computeRho(fu.linRangeByStep(10, 30))
    dcca.saveObject(get_object_path('dcca_obj'))
    n_rho_load = fu.getObjectMember(get_object_path('dcca_obj', ext=True), 'nRho')
    rho_load = fu.getObjectMember(get_object_path('dcca_obj', ext=True), 'rho')
    n_thr_load = fu.getObjectMember(get_object_path('dcca_obj', ext=True), 'nThr')
    conf_up_load = fu.getObjectMember(get_object_path('dcca_obj', ext=True), 'confUp')
    conf_down_load = fu.getObjectMember(get_object_path('dcca_obj', ext=True), 'confDown')
    assert np.array_equal(n_load, n)
    assert np.array_equal(F_load, F)
    assert np.array_equal(n_rho_load, n_rho)
    assert np.array_equal(rho_load, rho)
    assert np.array_equal(n_thr_load, [])
    assert np.array_equal(conf_up_load, [])
    assert np.array_equal(conf_down_load, [])
    # thresholds
    n_thr, conf_up, conf_down = dcca.rhoThresholds(len(ts), fu.linRangeByStep(10, 30), 10, 0.95)
    dcca.saveObject(get_object_path('dcca_obj'))
    n_thr_load = fu.getObjectMember(get_object_path('dcca_obj', ext=True), 'nThr')
    conf_up_load = fu.getObjectMember(get_object_path('dcca_obj', ext=True), 'confUp')
    conf_down_load = fu.getObjectMember(get_object_path('dcca_obj', ext=True), 'confDown')
    assert np.array_equal(n_load, n)
    assert np.array_equal(F_load, F)
    assert np.array_equal(n_rho_load, n_rho)
    assert np.array_equal(rho_load, rho)
    assert np.array_equal(n_thr_load, n_thr)
    assert np.array_equal(conf_up_load, conf_up)
    assert np.array_equal(conf_down_load, conf_down)


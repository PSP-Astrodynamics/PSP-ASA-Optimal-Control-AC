function[w_err] = w_error(test_vec, true_vec, test_w, true_w)


[~,quat_err] = quat_orien_err(test_vec,true_vec); 


w_err = true_w - quat_err % rotate the test w by quat error
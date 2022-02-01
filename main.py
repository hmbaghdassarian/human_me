from human_me.preprocess.correct_inputs import correct_model, correct_psim
from human_me.build.build_me_model import build_me

cm_1, cm_2, me_input_model = correct_model(model_file= ) 
corrected_psim_me, non_machinery, revised_genes = correct_psim(me_input_model, 
                                                                psim_df=build_local_path + 'psim_me.h5', # try with default and explicitly stating
                                                                non_machinery=None) 
me_model, me_builder = build_me(me_input_model=me_input_model, psim_me=corrected_psim_me, **kwargs)
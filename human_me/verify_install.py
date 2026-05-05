"""Verify successful human_me installation"""

import sys

def main():
    print("Testing imports")
    import human_me
    from human_me.data.file_paths import build_files_url, build_local_path, input_local_path
    from human_me.build.build_me_model import build_me
    from human_me.me_solver.solve_me import solve_lp
    from human_me.utils import downstream_analyses

    import cobra
    import cobra.test
    import numpy as np
    import pandas as pd

    print("Testing solver...")
    model = cobra.test.create_test_model("textbook")
    m_sln_q, stat, _ = solve_lp(
        me_model=model,
        mu_val=np.nan,
        objective={'Biomass_Ecoli_core': 1},
        verbosity=False,
        precision='quad',
    )

    assert stat.item() == 0, 'qMINOS solver cannot identify a solution'

    formatted = pd.DataFrame({'reaction_id': [r.id for r in model.reactions]})
    formatted['flux'] = formatted.reaction_id.apply(
        lambda r_id: m_sln_q[model.reactions.index(r_id)]
    )
    formatted.set_index('reaction_id', inplace=True)

    assert np.isclose(
        formatted.loc['Biomass_Ecoli_core', 'flux'],
        model.slim_optimize()
    ), 'qMINOS solver does not identify correct solution'

    print("✓ Installation verified successfully")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"✗ Verification failed: {e}", file=sys.stderr)
        sys.exit(1)
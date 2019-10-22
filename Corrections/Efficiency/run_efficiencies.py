#!/usr/bin/env python3

import click
import datetime
from subprocess import call
from functools import partial

shell_call = partial(call, shell=True)

def get_systematic_type(ctx, params, value):
    if value in ['Nominal', 'Binned']:
        return value
    elif ctx.params['collision_type'] == 'pPb':
        up_or_down = click.prompt("Upward or downward fluctuation?", default='Up')
        return value, up_or_down
    else:
        systematic_type = click.prompt("What kind of TnP to use for pp systematics?", default='Trigger')
        up_or_down = click.prompt("Upward or downward fluctuation?", default='Up')
        return value, up_or_down, systematic_type


@click.command()
@click.option('--identifier',
              default=str(datetime.date.today()).replace('-', '_'),
              prompt='Folder name',
              help='Call the folder something logical')
@click.option('--state', default='1S', prompt='Which Upsilon state?', help='Can be 1S, 2S, 3S or all_states')
@click.option('--collision_type', default='pPb', prompt='Collision type',
              help='Can be pp or pPb')      # Generalization to other collision types
@click.option('--tnp_correction', default='Nominal', prompt='What TnP scale factors to use',
              callback=get_systematic_type,
              help='Can be Nominal, Systematics, Binned (Binned is only for pp)')
def run_efficiencies(identifier, state, collision_type, tnp_correction):
    if (not isinstance(tnp_correction, tuple)):
        up_or_down = None
        systematic_type = None
    elif len(tnp_correction) == 2:
        tnp_correction, up_or_down = tnp_correction
        systematic_type = None
    elif len(tnp_correction) == 3:
        tnp_correction, up_or_down, systematic_type = tnp_correction
    else:
        raise ValueError('You did not enter the write options for this collision type')

    if up_or_down is not None:
        identifier += up_or_down

    if (collision_type == 'pp') and (systematic_type is not None):
        identifier += systematic_type

    shell_call(f'mkdir eff_{collision_type}{identifier}')
    if state in ['1S', '2S', '3S']:
        state = int(state[0])
        run_one_state(state, collision_type, tnp_correction, identifier, systematic_type, up_or_down)
    elif state == 'all_states':
        for tmp_state in [1, 2, 3]:
            run_one_state(tmp_state, collision_type, tnp_correction,
                          identifier, systematic_type, up_or_down)
    else:
        raise ValueError('Give it a state')


def run_one_state(state, collision_type, tnp_correction, identifier, systematic_type, up_or_down):
    if collision_type == 'pPb':
        ispPb = 1
    elif collision_type == 'pp':
        ispPb = 0
        if (tnp_correction == 'Systematics') and (systematic_type is None):
            raise ValueError(
                'You want systematic deviation due to tnp corrections for pp but did not specify the tnp type')
        if (tnp_correction == 'Systematics') and (up_or_down is None):
            raise ValueError('You forgot to specify upward or downward fluctuation')
    else:
        raise ValueError('You didnt specify pp or pPb')

    shell_call(f'cp dimuEff_RpA_ana.C dimuEff_oniaMode{state}_{collision_type}_{tnp_correction}_{identifier}.C')
    shell_call(f'sed -i "" "s/VVV/{state}/g;" dimuEff_oniaMode{state}_{collision_type}_{tnp_correction}_{identifier}.C')
    shell_call(f'sed -i "" "s/WWW/{ispPb}/g;" dimuEff_oniaMode{state}_{collision_type}_{tnp_correction}_{identifier}.C')
    shell_call(f'sed -i "" "s/XXX/{collision_type}/g;" dimuEff_oniaMode{state}_{collision_type}_{tnp_correction}_{identifier}.C')
    shell_call(f'sed -i "" "s/tnp_choice/{tnp_correction}/g;" dimuEff_oniaMode{state}_{collision_type}_{tnp_correction}_{identifier}.C')
    shell_call(f'sed -i "" "s/pp_tnp_sys_choice/{systematic_type or -1}/g;" dimuEff_oniaMode{state}_{collision_type}_{tnp_correction}_{identifier}.C')
    shell_call(f'sed -i "" "s/RpA_ana/oniaMode{state}_{collision_type}_{tnp_correction}_{identifier}/g;" dimuEff_oniaMode{state}_{collision_type}_{tnp_correction}_{identifier}.C')
    shell_call(f'sed -i "" "s/TAG/{identifier}/g;" dimuEff_oniaMode{state}_{collision_type}_{tnp_correction}_{identifier}.C')
    shell_call(f'sed -i "" "s/up_or_down/{up_or_down or -1}/g;" dimuEff_oniaMode{state}_{collision_type}_{tnp_correction}_{identifier}.C')
    shell_call(f'cp dimuEff_oniaMode{state}_{collision_type}_{tnp_correction}_{identifier}.C eff_{collision_type}{identifier}/dimuEff_oniaMode{state}_{collision_type}_{tnp_correction}_{identifier}.C')
    shell_call(f'root -l -q -b ./dimuEff_oniaMode{state}_{collision_type}_{tnp_correction}_{identifier}.C')
    shell_call(f'rm dimuEff_oniaMode{state}_{collision_type}_{tnp_correction}_{identifier}.C')


if __name__ == '__main__':
    run_efficiencies()

#!/usr/bin/env bash
[ -d 'SIMULATION_FILES' ] && rm -R SIMULATION_FILES 
[ -d 'RESTART' ]          && rm -R RESTART
[ -f 'OUTPUT' ]           && rm OUTPUT
[ -f 'SET_SIMULATION' ]   && rm SET_SIMULATION
[ -f 'SET_KPOINTS'    ]   && rm SET_KPOINTS
[ -f 'RECORD_MODELS'  ]   && rm RECORD_MODELS

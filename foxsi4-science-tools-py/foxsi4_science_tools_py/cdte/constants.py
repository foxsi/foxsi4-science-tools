"""Module to store constant information for CdTe."""

STRIPS_PER_SIDE = 128 # numbered 0-127
UNPHYSICAL_STRIP = 128 # this strip index doesn't match a physical strip
UNPHYSICAL_ENERGY = -1.0 # a value to use in energy arrays where there is not a valid value
UNPHYSICAL_ENERGY_DIFF = -100.0 # a value to use in energy difference arrays where there is not a valid value
UNPHYSICAL_NHIT = -1 # a dummy value for a trigger entry that doesn't have a physical number of hits
UNPHYSICAL_SPLINE = -1.0e+6 # a dummy value for a spline energy value
# mappings of different strip number definitions to their ASICs, ranges are inclusive
ASICS_TO_STRIPS_MAPPING = {"INDEX_PT":{"ASIC0":[0, 63], "ASIC1":[64, 127]},
                           "INDEX_AL":{"ASIC3":[0, 63], "ASIC4":[64, 127]},
                           "INDEX_AL_ALTERNATIVE":{"ASIC3":[128, 191], "ASIC4":[192, 255]},
                           "REMAPCH_PT":{"ASIC0":[0, 63], "ASIC1":[64, 127]},
                           "REMAPCH_AL":{"ASIC3":[64, 127], "ASIC4":[0, 63]},
                           "REMAPCH_AL_ALTERNATIVE":{"ASIC3":[192, 255], "ASIC4":[128, 191]},
                           }
# mapping the strip re-mapped values to their pitch regions, ranges are inclusive
STRIP_TO_PITCH_MAPPING = {"GUARD_RING_A":[0, 3], # guard ring, 100 um pitch
                          "REGION2_A":[4, 27], # 100 um pitch
                          "REGION1_A":[28, 47], # 80 um pitch
                          "REGION0":[48, 78], # 60 um pitch
                          "REGION1_B":[79, 98], # 80 um pitch
                          "REGION2_B":[99, 122], # 100 um pitch
                          "GUARD_RING_B":[123, 127], # guard ring, 100 um pitch
                          }

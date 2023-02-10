hbll_highlights <- c(
  # "petrale sole",
  "arrowtooth flounder",
  "big skate",
  "canary rockfish",
  "china rockfish", # 13% positive
  "copper rockfish", # 11% positive
  "lingcod",
  "longnose skate",
  "north pacific spiny dogfish",
  "pacific cod",
  "quillback rockfish",
  "redbanded rockfish",
  "rosethorn rockfish",
  "rougheye/blackspotted rockfish complex",
  "sandpaper skate", # 7.6 % positive, gets filtered out
  "shortspine thornyhead", # only 6.7% of samples are positive, 8% pos restricted
  "silvergray rockfish",
  "southern rock sole",
  "spotted ratfish",
  "tiger rockfish",
  "yelloweye rockfish"
) %>%
  tolower() %>%
  sort() %>%
  unique()

syn_highlights <- c(
  "Harlequin Rockfish", ####
  "Pacific Sanddab", #### abundant flatfish
  "Sand Sole", #### abundant flatfish
  "Big Skate",
  "Longnose Skate",
  "sandpaper skate", ###
  "North Pacific Spiny Dogfish",
  "Spotted Ratfish",
  "kelp greenling", ###
  "Lingcod",
  "Pacific Cod",
  "Walleye Pollock",
  "Bocaccio",
  "Canary Rockfish",
  "darkblotched rockfish", ###
  "greenstriped rockfish", ###
  "Pacific Ocean Perch",
  "Redstripe Rockfish",
  "Redbanded Rockfish",
  "rougheye/blackspotted rockfish complex",
  "sharpchin rockfish", ###
  "Shortraker Rockfish", # not fitting well! almost no data when restricted (9% pos)
  "Shortspine Thornyhead",
  "Silvergray Rockfish",
  "Widow Rockfish",
  "Yellowmouth Rockfish",
  "Yellowtail Rockfish",
  "Arrowtooth Flounder",
  "Butter Sole", ###
  "Curlfin Sole",
  "Dover Sole",
  "English Sole",
  "Flathead Sole",
  "Petrale Sole",
  "Rex Sole",
  "Slender Sole",
  "southern rock sole"
) %>%
  tolower() %>%
  sort() %>%
  unique()

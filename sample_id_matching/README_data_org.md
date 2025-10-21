# HGSOC data organization

## Original data organization
Two folders in `/pl/active/cgreene-sc-hgsoc/`:
```
penn_HGSOC/
├── basespace_md5sum.txt
├── HGSOC-2216_ds.88ec815863be40ad8260ac629bb698ad
├── HGSOC-2221_ds.50a404bc0a584d6ba76bbf39fd422af4
├── HGSOC-2230_ds.89b2b7d602bc4689843d1d90fcca4a5f
├── HGSOC-2238_ds.1e19d229b15747be9f77e064ebe8a467
├── HGSOC-2246_ds.782ffae2e3b84349b92a0914dbd1c60f
├── HGSOC-2249_ds.67913644ef28455eabd1e6a017f2bd0b
├── HGSOC-2268_ds.083794fba95c4554b7d41ccc652ac943
├── HGSOC-2278_ds.e70b101c0d594d709e7ebf0367726bb4
├── HGSOC-2296_ds.9ae05fc537eb44f8bd1670dfa2df9688
├── HGSOC-2364_ds.68ca287afa0444749d30334394b5baf4
├── HGSOC-2408_ds.41a1186148774fdb956a9308ee7060f9
├── HGSOC-2416_ds.6d71d835e53a45b6968809911623588e
├── HGSOC-2423_ds.31453df6b8844dd091adb7e279f9e30f
├── HGSOC-2430_ds.97f1621f80a041d2975059fdff96ca4d
├── HGSOC-2455_ds.39c0e8300e5b497f8396013a428bc05d
├── HGSOC-2460_ds.e975eb3f47824d4ab8c8b30e6feca083
├── HGSOC-2466_ds.8a9f697f122c48a89afb18fe7f7cb16f
├── HGSOC-2477_ds.97607d06c5cd43cc8381dd5f742d6e77
├── HGSOC-2507_ds.22c8261caaf246b0beb0f8520390fb3f
├── HGSOC-2514_ds.5046358a6493466790c134defe660ce0
├── HippenA_230626_A00901_0948_BHFGW3DRX3_392180799.json
└── md5sum.txt
```
and
```
ariel_sc_HGSOC/
├── GreeneC_OrmanD_scV3_HM5CCDSX5_230201
│   ├── fastq
│   ├── Pool10-CSP-CaseyGreene-01052023
│   ├── Pool10-GEX-CaseyGreene-01052023
│   ├── Pool1-CSP-CaseyGreene-11302022
│   ├── Pool1-GEX-CaseyGreene-11302022
│   ├── Pool2-CSP-CaseyGreene-11302022
│   ├── Pool2-GEX-CaseyGreene-11302022
│   ├── Pool3-CSP-CaseyGreene-12122022
│   ├── Pool3-GEX-CaseyGreene-12122022
│   ├── Pool4-CSP-CaseyGreene-12122022
│   ├── Pool4-GEX-CaseyGreene-12122022
│   ├── Pool5-CSP-CaseyGreene-01042023
│   ├── Pool5-GEX-CaseyGreene-01042023
│   ├── Pool6-CSP-CaseyGreene-01042023
│   ├── Pool6-GEX-CaseyGreene-01042023
│   ├── Pool7-CSP-CaseyGreene-01052023
│   ├── Pool7-GEX-CaseyGreene-01052023
│   ├── Pool8-CSP-CaseyGreene-01052023
│   ├── Pool8-GEX-CaseyGreene-01052023
│   ├── Pool9-CSP-CaseyGreene-01052023
│   └── Pool9-GEX-CaseyGreene-01052023
├── HippenA_230414_A00901_0909_BH3FJ2DRX3
│   ├── HGSOC-2018_ds.9c223da745e04293a5565413e51ee410
│   ├── HGSOC-2023_ds.44924ac21313411bb3955a62cbcc162a
│   ├── HGSOC-2126_ds.793e5acb10f1443dbc84d0c67799d054
│   ├── HGSOC-2129_ds.d31c91669e554f5c8f2a9158dfc5a932
│   ├── HGSOC-2202_ds.bc51053465b24dd29dd404f24b18bc5b
│   ├── HGSOC-2209_ds.42d44c8ef289402e88c547b92921b065
│   ├── HGSOC-2221_ds.887f415a92e44961bff818d3a906477b
│   ├── HGSOC-2238_ds.e7767fb943124326a98f3affeb9d722c
│   ├── HGSOC-2240_ds.4ca7d3c1d7f944bebf319f3504a7ebe1
│   ├── HGSOC-2249_ds.a255da8cd9d845b69a3773c923d22a21
│   ├── HGSOC-2278_ds.7d4db1da3e5c4483988e41030b8f81ad
│   ├── HGSOC-2313_ds.a801e1416260459aaf2a1c3c3a5f5133
│   ├── HGSOC-2364_ds.4e10fcd4ca1144a7a18bb1aa7cff4aab
│   ├── HGSOC-2401_ds.d913e45f78bf4299a16d2b713937961d
│   ├── HGSOC-2407_ds.b187a51e716a4100ac22e5fb8a2b9625
│   ├── HGSOC-2408_ds.ef75a27f504940d2b2390241b0e3abe3
│   ├── HGSOC-2416_ds.1dcfade89cc644baad1faaac5c89573a
│   ├── HGSOC-2460_ds.8aea59922ff64a0d9e7da88b1593b72c
│   ├── HGSOC-2477_ds.0af6a47fc92a4675b1f6c60ae2a863d9
│   ├── HGSOC-2483_ds.2007dc0d3e9048f0b64e4792f533e622
│   ├── HGSOC-2507_ds.35900512328043c7ad14c0b7faa6514b
│   ├── HGSOC-2526_ds.22414930847544d28aeb703fe56f43cf
│   ├── HippenA_230414_A00901_0909_BH3FJ2DRX3_386708325.json
│   └── README
├── HippenA_230418_A00901_0912_BHK72LDMXY
│   ├── HGSOC-2094_ds.a994179dcb21402e8e8b71510093e568
│   ├── HGSOC-2186_ds.18debb70861d4f0b8cbe8f15664d40e6
│   ├── HGSOC-2216_ds.444bcb3644fd43b78470e7a175e23a53
│   ├── HGSOC-2230_ds.dc1081facbaa4199b207fd619422c83b
│   ├── HGSOC-2246_ds.42d245f9cdda411fba8c77fb78581984
│   ├── HGSOC-2309_ds.9942ee68daea4ab3ab85bfdb371b0029
│   ├── HGSOC-2423_ds.bfb3e1cf3f8b4448820689769b0ce937
│   ├── HGSOC-2455_ds.17f92227aa344ad7851c0da5c10baddf
│   ├── HGSOC-2466_ds.08b3c08e64a44ef4a48c35baa4796d1b
│   ├── HGSOC-2514_ds.2c962c8e2630472692b1e86d4d389b43
│   ├── HippenA_230418_A00901_0912_BHK72LDMXY_386762381.json
│   └── README
├── HippenA_230509_A00901_0926_AHK3THDMXY
│   ├── HGSOC-2018_ds.6468f190988c4f3cb7666384e5d5c43b
│   ├── HGSOC-2023_ds.77dbaa79198c4ea29103d18def53d880
│   ├── HGSOC-2094_ds.b549d17d0c95422caa89063655dc7595
│   ├── HGSOC-2126_ds.5f084153ce074da880f2fa5661c000c7
│   ├── HGSOC-2129_ds.eebdcd4b684b488ba880d33dab05a78c
│   ├── HGSOC-2186_ds.b07ea1855e0847e0b289f8c733eaadd7
│   ├── HGSOC-2202_ds.44cc7dc92f4f41d6b19aad9120edf9e7
│   ├── HGSOC-2209_ds.d4f2f480ab1144a09ead83b5ff3fd199
│   ├── HGSOC-2240_ds.b9f35f116b2e4247bd6bfdb9841663ad
│   ├── HGSOC-2309_ds.20b83c26912b4789adf0ba742df749c7
│   ├── HGSOC-2313_ds.ab79894e0bda40089d6546a599794b4f
│   ├── HGSOC-2401_ds.f9c04c189deb4200933e4e873974f5c3
│   ├── HGSOC-2407_ds.92eab3b5005e4209b8323e903cd73e27
│   ├── HGSOC-2444_ds.fdd3a47b8c5a445f9918dd57215051f3
│   ├── HGSOC-2483_ds.ad6f8e7115f94586a9b86d217829ee76
│   ├── HGSOC-2526_ds.d8d90e02e98d4986bb06fb5256a51e13
│   └── HippenA_230509_A00901_0926_AHK3THDMXY_388260874.json
└── md5sums.txt
```

## What's where in the original organization scheme?

- Pooled single-cell RNA-seq: 
    - `/ariel_sc_HGSOC/GreeneC_OrmanD_scV3_HM5CCDSX5_230201`
    - fastq files in the subfolder `fastq`.
    - various CellRanger outputs in the subfolders labelled with pool numbers, e.g. `Pool1-GEX-CaseyGreene-11302022`
    - a total of 37 samples
- Bulk dissociated poly A+ RNA-seq: 
    - `/ariel_sc_HGSOC/HippenA_230414_A00901_0909_BH3FJ2DRX3` and 
    - `/ariel_sc_HGSOC/HippenA_230418_A00901_0912_BHK72LDMXY`
    - fastq files in subfolders labelled with sample ids
    - a total of 32 samples
- Bulk chunk ribo (aka rRNA-) RNA-seq: 
    - `/ariel_sc_HGSOC/HippenA_230509_A00901_0926_AHK3THDMXY` and 
    - `/penn_HGSOC`
    - fastq files in subfolders labelled with sample ids
    - a total of 36 samples

Code and notes from before we figured how the data was organized often refers to the data by dates
- 230414 (or 0414) and 230418 or (0418) for the bulk dissociated poly A+ data
- 230509 (or 0509) and 230626 or (0626) for the bulk chunk ribo data

### How do we know this?

A combination of old README's and differential expression analysis, see [hgsoc_data_detective](https://github.com/greenelab/hgsoc_data_detective) (Github repo).

### Sample (mis)match note

In principle, the different RNA-seq modalities listed above have been executed on samples from the same donors and are labelled as such. 
However, a significant number of samples have been mislabelled with respect to each other. 

For now, the samples remain labelled with their original (sometimes wrong) sample ids. This is being addressed in the [hgsoc_power](https://github.com/greenelab/hgsoc-power) project (sub project / folder `sample_id_matching`).


## Re-organization (Oct 2025)

I would like to reorganize the data in the two folders described above. It will all stay in `/pl/active/cgreene-sc-hgsoc/`, but the organization within that would change:

```
ariel_sc_HGSOC/
├── pooled_sc/
|   ├── fastq/
|   ├── poolX-subfolders/
├── bulk_diss_polyA/
├── bulk_chunk_ribo/
```

The main changes are combining a few folders together and renaming them. 
Also, I want to take out the layer of folders that have all the sample ids in them, since the sample ids are also in the fastq file names themselves.

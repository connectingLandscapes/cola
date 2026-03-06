import ee



def l4to9_c2_select_rename(img):
    """

    :param img: landsat 4 to 9 collection 2 image
    :return: landsat image with renamed bands
    """

    l457_bandnames = ee.List(['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7'])
    l89_bandnames = ee.List(['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7'])

    bandname_dict = ee.Dictionary({
        'LANDSAT_4': l457_bandnames,
        'LANDSAT_5': l457_bandnames,
        'LANDSAT_7': l457_bandnames,
        'LANDSAT_8': l89_bandnames,
        'LANDSAT_9': l89_bandnames})

    band_names = ee.List(bandname_dict.get(img.get('SPACECRAFT_ID'))).cat(ee.List(['QA_PIXEL', 'QA_RADSAT']))
    new_band_names = ee.List(['blue', 'green', 'red', 'NIR', 'SWIR1', 'SWIR2', 'QA_PIXEL', 'QA_RADSAT'])
    img_select_rename = img.select(band_names, new_band_names)

    return img_select_rename


def l4to9_c2_scaleoff(img):
    """

    :param img: landsat 4 to 9 collection 2 image
    :return: landsat image scaled to units of surface reflectance [0-1]
    """

    sel_rn_img = l4to9_c2_select_rename(img)
    sr_img = (sel_rn_img.select(['blue', 'green', 'red', 'NIR', 'SWIR1', 'SWIR2'])
                        .multiply(0.0000275)
                        .add(-0.2)
                        .clamp(0,1)
                        .addBands(img.select(['QA_PIXEL', 'QA_RADSAT']))
                        .copyProperties(img, img.propertyNames()))

    return sr_img


def l4to9_c2_qa_mask_clouds(img):

    def extract_qa_bits(qa_band, bit_start, bit_end):
        """

        :param qa_band: landsat 4 to 9 collection 2 image PIXEL_QA band
        :param bit_start: starting bit in qa_band
        :param bit_end: ending bit in qa_band
        :return: a transformed quality band
        """

        num_bits = bit_end - bit_start + 1
        qa_bits = qa_band.rightShift(bit_start).mod(ee.Number(2).pow(num_bits))
        return qa_bits

    # extract relevant parts of the quality mask
    qa = img.select(['QA_PIXEL'])
    cloud_mask = qa.bitwiseAnd(1 << 3).eq(0)
    cirrus_mask = extract_qa_bits(qa, 14, 15).lte(2)
    dilatedcloud_mask = qa.bitwiseAnd(1 << 1).eq(0)
    cloudshadow_mask = qa.bitwiseAnd(1 << 4).eq(0)

    # radiometric saturation mask
    saturation_mask = img.select('QA_RADSAT').eq(0)

    combined_mask = cloud_mask.And(cloudshadow_mask).And(cirrus_mask).And(saturation_mask)

    return img.updateMask(combined_mask)


def l4to9_c2_indices(img):
    """

    :param img: landsat 4 to 9 collection 2 image
    :return: landsat image with spectral indices added
    """

    # Normalized Difference Vegetation Index (NDVI)
    # Ref: Tucker, 1979
    NDVI = img.normalizedDifference(['NIR', 'red']).rename('NDVI')

    # Normalized Difference Moisture Index (NDMI)
    NDMI = img.normalizedDifference(['NIR', 'SWIR1']).rename('NDMI')

    # Normalized Burn Ratio
    NBR = img.normalizedDifference(['NIR', 'SWIR2']).rename('NBR')

    # Normalized Burn Ratio 2
    NBR2 = img.normalizedDifference(['SWIR1', 'SWIR2']).rename('NBR2')

    return ee.Image.cat(img, NDVI, NDMI, NBR, NBR2)

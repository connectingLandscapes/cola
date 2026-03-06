import ee
import math


def mosaic_ccdc_img_tiles(ccdc_result_path:str = None):
    """Mosaic CCDC array images

    :param ccdc_result_path: GEE asset path with CCDC results
    :return: a mosaiced array image of CCDC results
    """

    ccdc_asset_list = [i['id'] for i in ee.data.listAssets(ccdc_result_path)['assets']]
    ccdc_img = ee.ImageCollection(ccdc_asset_list).mosaic()

    return ccdc_img


def normalize_intercept_band(ccdc_img: ee.Image = None,
                             band: str = None):
    """Normalize the y-intercept of a band relative to the middle date

    :param ccdc_img: CCDC array image
    :param band: band name
    :return: a CCDC image with normalization applied to a band
    """

    mid = ccdc_img.select('tStart').add(ccdc_img.select('tEnd')).divide(2).toArray(1)
    coefs = ccdc_img.select(ee.String(band).cat(ee.String('_coefs')))
    adjusted = coefs.arraySlice(1, 0, 1).add(coefs.arraySlice(1, 1, 2).multiply(mid))

    return adjusted.arrayCat(coefs.arraySlice(1, 1), 1).rename(ee.String(band).cat(ee.String('_coefs')))


def normalize_intercept(ccdc_img: ee.Image = None,
                        bands: list = None):
    """Normalize the y-intercept of bands relative to the middle date

    :param ccdc_img: CCDC array image
    :param bands: list of band names
    :return: a CCDC image with normalization applied to select bands
    """

    result = []
    for band in bands:
        result.append(normalize_intercept_band(ccdc_img=ccdc_img,
                                               band=band))

    return ccdc_img.addBands(ee.Image(result), None, True)


# TODO: possible to get closest segment in addition to previous or next?
def get_segment_coefs_band(ccdc_img: ee.Image = None,
                           year_frac: float = None,
                           band: str = None,
                           use_next: bool = False):
    """Get CCDC segment coefficients for a band

    :param ccdc_img: CCDC array image
    :param year_frac: decimal year
    :param band: band name
    :param use_next: if a date is between segments, use next segment (True) or previous (False)
    :return: ee.Image of CCDC segment coefficients for a band
    """

    tStart = ccdc_img.select('tStart')
    tEnd = ccdc_img.select('tEnd')
    coefs = ccdc_img.select(ee.String(band).cat(ee.String('_coefs')))

    # use previous segments if requested year_frac is not belong to a segment.
    coefNames = ["INTP", "SLP", "COS", "SIN", "COS2", "SIN2", "COS3", "SIN3"]
    coefNames = ['ccsegco_' + band + '_' + i.lower() for i in coefNames]

    # use next segment
    vicNext = tEnd.gt(year_frac)
    #//is year_frac pass the end of time series
    fixNext = vicNext.arrayReduce(ee.Reducer.sum(), [0]).arrayFlatten([['count']])
    selectedNext = coefs.arrayMask(vicNext.toArray(1)).arraySlice(0, 0, 1).arrayProject([1]).arrayFlatten([coefNames])

    # use previous segment
    vicPrev = tStart.lte(year_frac)
    #//is year_frac before the beginning of time series
    fixPrev = vicPrev.arrayReduce(ee.Reducer.sum(), [0]).arrayFlatten([['count']])
    selectedPrev = coefs.arrayMask(vicPrev.toArray(1)).arraySlice(0, -1).arrayProject([1]).arrayFlatten([coefNames])

    # fix beginning and end of time series issue
    fixedSelectedNext = selectedNext.where(fixNext.eq(0), selectedPrev)
    fixedSelectedPrev = selectedPrev.where(fixPrev.eq(0), selectedNext)
    if use_next:
      return fixedSelectedNext
    else:
      return fixedSelectedPrev


def get_segment_coefs(ccdc_img: ee.Image = None,
                      year_frac: float = None,
                      bands: list = None,
                      normalize: bool = True,
                      use_next: bool = False):
    """Get CCDC segment coefficients for multiple bands

    :param ccdc_img: CCDC array image
    :param year_frac: decimal year
    :param bands: list of band names
    :param normalize: normalize the y-intercept of bands relative to the middle date
    :param use_next: if a date is between segments, use next segment (True) or previous (False)
    :return: ee.Image of CCDC segment coefficients for multiple bands
    """

    if normalize:
        ccdc_img = normalize_intercept(ccdc_img=ccdc_img,
                                       bands=bands)

    result = []
    for band in bands:
        result.append(get_segment_coefs_band(ccdc_img=ccdc_img,
                                             year_frac=year_frac,
                                             band=band,
                                             use_next=use_next))

    return ee.Image(result).float()


# TODO: possible to extrapolate a segment? KO looking into Wiell
def get_synthetic_band(ccdc_img: ee.Image = None,
                       year_frac: float = None,
                       band: str = None,
                       use_next: bool = False):
    """Get CCDC synthetic image for a band

    :param ccdc_img: CCDC array image
    :param year_frac: decimal year
    :param band: band name
    :param use_next: if a date is between segments, use next segment (True) or previous (False)
    :return: ee.Image of CCDC synthetic (fitted) values for a band
    """

    omega = 2.0 * math.pi
    year_frac = ee.Number(year_frac)
    imageT = ee.Image.constant([1, year_frac,
                                year_frac.multiply(omega).cos(),
                                year_frac.multiply(omega).sin(),
                                year_frac.multiply(omega * 2).cos(),
                                year_frac.multiply(omega * 2).sin(),
                                year_frac.multiply(omega * 3).cos(),
                                year_frac.multiply(omega * 3).sin()])

    coefs = get_segment_coefs_band(ccdc_img=ccdc_img,
                                   year_frac=year_frac,
                                   band=band,
                                   use_next=use_next)

    new_band_name = ee.String('ccsynth_').cat(ee.String(band))

    return imageT.multiply(coefs).reduce(ee.Reducer.sum()).rename(new_band_name)


def get_synthetic(ccdc_img: ee.Image = None,
                  year_frac: float = None,
                  bands: list = None,
                  use_next: bool = False):
    """Get CCDC synthetic image for multiple bands

    :param ccdc_img: CCDC array image
    :param year_frac: decimal year
    :param bands: list of band names
    :param use_next: if a date is between segments, use next segment (True) or previous (False)
    :return: ee.Image of CCDC synthetic (fitted) values for multiple bands
    """

    result = []
    for band in bands:
        result.append(get_synthetic_band(ccdc_img=ccdc_img,
                                         year_frac=year_frac,
                                         band=band,
                                         use_next=use_next))

    return ee.Image(result)


def get_segment_phase_amp(coefs_img: ee.Image = None,
                          sinExpr: str = '.*sin.*',
                          cosExpr: str = '.*cos.*'):
    """Get phase and amplitude of CCDC segment

    :param coefs_img: CCDC segment coefficient image
    :param sinExpr: regular expression to use for getting SIN coefficients
    :param cosExpr: regular expression to use for getting COS coefficients
    :return: ee.Image with phase, amplitude, and amplitude/phase ratio
    """

    sin = coefs_img.select(sinExpr)
    cos = coefs_img.select(cosExpr)

    #// Scale to [0, 1] from radians.
    phase = sin.atan2(cos) \
      .unitScale(-3.14159265359, 3.14159265359) \
      .multiply(365) #// To get phase in days!

    amplitude = sin.hypot(cos)
    
    ap_ratio = amplitude.divide(phase)

    apNames = ap_ratio.bandNames().map(lambda x: ee.String(x).replace('_sin', '_amphratio'))
    phaseNames = phase.bandNames().map(lambda x: ee.String(x).replace('_sin', '_phase'))
    amplitudeNames = amplitude.bandNames().map(lambda x: ee.String(x).replace('_sin', '_amp'))

    return phase.rename(phaseNames).addBands(amplitude.rename(amplitudeNames)).addBands(ap_ratio.rename(apNames))


def get_largest_mag_tbreak_thresh_img(ccdc_img: ee.Image = None,
                                      start_year: int = None,
                                      end_year: int = None,
                                      band: ee.String = None,
                                      mag_thresh_abs: ee.Number = None):
    """Get the year of largest CCDC segment break greater than an absolute magnitude threshold

    :param ccdc_img: CCDC array image
    :param start_year: starting year for counting segment breaks
    :param end_year: ending year for counting segment breaks
    :param band: band name to use for determining a break
    :param mag_thresh_abs: absolute value of change to use for determining a break
    :return: the year of absolute largest break for a band
    """

    band_mag_name = ee.String(band).cat(ee.String('_magnitude'))
    band_mag = ccdc_img.select(band_mag_name)
    tBreak = ccdc_img.select('tBreak')
    mask = (tBreak.lt(ee.Number(end_year))
                  .And(tBreak.gte(ee.Number(start_year))
                  .And(band_mag.abs().gte(mag_thresh_abs))).toArray(0))
    dates = tBreak.arrayMask(mask).arrayPad([1])
    magnitudes = ccdc_img.select(band_mag_name).arrayMask(mask).arrayPad([1])
    maxIndex = magnitudes.abs().arrayArgmax().arrayFlatten([['index']])
    sel_mag = magnitudes.arrayGet(maxIndex)
    sel_tBreak = dates.arrayGet(maxIndex).selfMask()

    return sel_tBreak


def get_time_since_largest_mag_break(ccdc_img:ee.Image = None,
                                     start_year:int = None,
                                     year_frac: float = None,
                                     band: ee.String = None,
                                     mag_thresh_abs: ee.Number = None):
    """Get the years since the largest CCDC segment break greater than an absolute magnitude threshold

    :param ccdc_img: CCDC array image
    :param start_year: starting year for counting segment breaks
    :param year_frac: fractional year
    :param band: band name to use for determining a break
    :param mag_thresh_abs: absolute value of change to use for determining a break
    :return: years since the year of absolute largest break for a band
    """

    # time since break above magnitude
    year_break = get_largest_mag_tbreak_thresh_img(ccdc_img=ccdc_img,
                                                  start_year=start_year,
                                                  end_year=year_frac,
                                                  band=band,
                                                  mag_thresh_abs=mag_thresh_abs)
    new_band_name = ee.String('ccbreak_tsince_').cat(ee.String(band).toLowerCase())
    time_since_break = ee.Image(year_frac).subtract(year_break.unmask(start_year)).rename(new_band_name)

    return time_since_break


def get_time_since_largest_mag_break_bands(ccdc_img: ee.Image = None,
                                           start_year: int = None,
                                           year_frac: float = None,
                                           bands_list: list = None,
                                           thresh_list: list = None):
    """

    :param ccdc_img: CCDC array image
    :param start_year: starting year for counting segment breaks
    :param year_frac: fractional year
    :param bands_list: list of band names to use for determining a break
    :param thresh_list: list of absolute values of change to use for determining a break, will be paired with band list
    :return:
    """

    result = []
    for i in range(0,len(bands_list)):
        band = bands_list[i]
        thresh = thresh_list[i]
        result.append(get_time_since_largest_mag_break(ccdc_img=ccdc_img,
                                                       start_year=start_year,
                                                       year_frac=year_frac,
                                                       band=band,
                                                       mag_thresh_abs=thresh))

    return ee.Image(result)


def get_n_breaks_img(ccdc_img: ee.Image = None,
                     start_year: int = None,
                     year_frac: int = None):
    """Get the number of CCDC segment breaks during a period

    :param ccdc_img: CCDC array image
    :param start_year: starting year for counting segment breaks
    :param year_frac: ending fractional year for counting segment breaks
    :return: ee.Image of the number of segments breaks during a period
    """

    tBreak = ccdc_img.select('tBreak')
    vic = tBreak.lt(ee.Number(year_frac)).And(tBreak.gte(ee.Number(start_year))).toArray(0)

    return tBreak.arrayMask(vic).arrayLength(0)
import loadErfa from './erfa.js';

const PI = Math.PI;
const DEG2RAD = PI / 180.0;
const RAD2DEG = 180.0 / PI;
const DAYSEC = 86400.0;

let ERFA = null;
let iersData = null;

// ERFA C-function wrappers
let eraJd2cal = null;
let eraDtf2d = null;
let eraAtoc13 = null;

/**
 * Loads IERS Bulletin A table for UT1-UTC and polar motion corrections.
 */
export async function loadIERSData() {
  if (!window.IERS_DATA) {
    throw new Error('[IERS Error] Could not find "window.IERS_DATA". Ensure iers.js is loaded in index.html.');
  }

  if (!window.IERS_DATA.data) {
    throw new Error('[IERS Error] Invalid format in "window.IERS_DATA". Missing "data" key.');
  }

  iersData = window.IERS_DATA.data;
  // return iersData;
}

/**
 * Returns interpolated IERS parameters given a Julian Date or Unix timestamp.
 * 
 * @param {number} timeInput - Either Unix Time (sec) or Julian Date (JD > 2000000)
 * @returns {{dUT1: number, xp: number, yp: number, isExact: boolean}}
 */
export function getIERS(timeInput) {
  if (!iersData) {
    console.warn("[IERS Warning] Table not loaded. Using zero corrections (dUT1=0, xp=0, yp=0).");
    return { dUT1: 0.0, xp: 0.0, yp: 0.0, isExact: false };
  }

  // Determine exact MJD whether input is Unix Time or Julian Date
  const exactMJD = timeInput > 2000000.0 
    ? timeInput - 2400000.5 
    : (timeInput / DAYSEC) + 40587.0;

  const mjd1 = Math.floor(exactMJD);
  const mjd2 = mjd1 + 1;
  const frac = exactMJD - mjd1;

  const rec1 = iersData[mjd1];
  const rec2 = iersData[mjd2];

  const extract = (rec, key) => {
    if (typeof rec === 'object' && rec !== null) return rec[key] ?? 0.0;
    if (key === 'dut1' && typeof rec === 'number') return rec;
    return 0.0;
  };

  if (rec1 && rec2) {
    return {
      dUT1: extract(rec1, 'dut1') + (extract(rec2, 'dut1') - extract(rec1, 'dut1')) * frac,
      xp:   extract(rec1, 'xp')   + (extract(rec2, 'xp')   - extract(rec1, 'xp'))   * frac,
      yp:   extract(rec1, 'yp')   + (extract(rec2, 'yp')   - extract(rec1, 'yp'))   * frac,
      isExact: true
    };
  }

  console.warn(`[IERS Warning] MJD ${mjd1} not found in Bulletin A table. Using uncorrected values.`);
  return { dUT1: 0.0, xp: 0.0, yp: 0.0, isExact: false };
}

/**
 * Initializes ERFA WebAssembly bindings and loads IERS Bulletin A table.
 */
export async function initErfa() {
  if (!ERFA) {
    ERFA = await loadErfa();
    eraJd2cal = ERFA.cwrap('eraJd2cal', 'number', ['number', 'number', 'number', 'number', 'number', 'number']);
    eraDtf2d  = ERFA.cwrap('eraDtf2d',  'number', ['string', 'number', 'number', 'number', 'number', 'number', 'number', 'number', 'number']);
    eraAtoc13 = ERFA.cwrap('eraAtoc13', 'number', [
      'string', 'number', 'number', 'number', 'number', 
      'number', 'number', 'number', 'number', 'number', 
      'number', 'number', 'number', 'number', 'number', 
      'number', 'number'
    ]);
  }

  // Load or refresh the IERS Bulletin A table
  await loadIERSData();
}

/**
 * Transforms topocentric horizontal coordinates to equatorial ICRS using ERFA.
 */
export function getEquatorialERFA(mode, jd, lat_deg, lon_deg, elevation, zen_deg, az_deg) {
  if (!ERFA) {
    throw new Error("[ERFA Error] Module not initialized. Call 'await initErfa()' first.");
  }

  // Wasm Memory Allocation
  const iyPtr   = ERFA._malloc(4);
  const imPtr   = ERFA._malloc(4);
  const idPtr   = ERFA._malloc(4);
  const fdPtr   = ERFA._malloc(8);
  const utc1Ptr = ERFA._malloc(8);
  const utc2Ptr = ERFA._malloc(8);
  const rcPtr   = ERFA._malloc(8);
  const dcPtr   = ERFA._malloc(8);

  try {
    // 1. Split Julian Date (iau_JD2CAL)
    const DJ1 = 2400000.5;
    const DJ2 = jd - 2400000.5;

    let status = eraJd2cal(DJ1, DJ2, iyPtr, imPtr, idPtr, fdPtr);
    if (status === -1) throw new Error("Error: unacceptable date in eraJd2cal");

    const iy = ERFA.getValue(iyPtr, 'i32');
    const im = ERFA.getValue(imPtr, 'i32');
    const id = ERFA.getValue(idPtr, 'i32');
    const fd = ERFA.getValue(fdPtr, 'double');

    const daySecTotal = fd * DAYSEC;
    const ihr = Math.floor(daySecTotal / 3600.0);
    const imn = Math.floor((daySecTotal % 3600.0) / 60.0);
    const sec = (daySecTotal % 3600.0) % 60.0;

    // 2. Compute 2-part UTC (iau_DTF2D)
    status = eraDtf2d('UTC', iy, im, id, ihr, imn, sec, utc1Ptr, utc2Ptr);
    if (status < 0) throw new Error(`Error: eraDtf2d failed with code ${status}`);

    const utc1 = ERFA.getValue(utc1Ptr, 'double');
    const utc2 = ERFA.getValue(utc2Ptr, 'double');

    // 3. Obtain Earth Orientation Parameters (EOP) for precise mode
    let dUT1 = 0.0;
    let xp = 0.0;
    let yp = 0.0;

    if (mode === 'P') {
      const eop = getIERS(jd);
      dUT1 = eop.dUT1;
      // Convert arcseconds (from Bulletin A) to radians
      xp = (eop.xp / 3600.0) * DEG2RAD;
      yp = (eop.yp / 3600.0) * DEG2RAD;
    }

    // 4. Coordinates & Atmospheric parameters
    const tp   = 'A';
    const ob1  = az_deg * DEG2RAD;
    const ob2  = zen_deg * DEG2RAD;
    const elong = lon_deg * DEG2RAD;
    const phi  = lat_deg * DEG2RAD;
    const hm   = elevation; // Height in meters above reference ellipsoid

    // 5. Transform (iau_ATOC13)
    status = eraAtoc13(
      tp, ob1, ob2, utc1, utc2, dUT1,
      elong, phi, hm, xp, yp, 0.0, 0.0, 0.0, 0.0,
      rcPtr, dcPtr
    );
    if (status < 0) throw new Error(`Error: eraAtoc13 failed with code ${status}`);

    const rc = ERFA.getValue(rcPtr, 'double');
    const dc = ERFA.getValue(dcPtr, 'double');

    let ra_deg = rc * RAD2DEG;
    if (ra_deg < 0) ra_deg += 360.0;
    const dec_deg = dc * RAD2DEG;

    return { ra_deg, dec_deg };

  } finally {
    // Memory release
    ERFA._free(iyPtr); ERFA._free(imPtr); ERFA._free(idPtr); ERFA._free(fdPtr);
    ERFA._free(utc1Ptr); ERFA._free(utc2Ptr); ERFA._free(rcPtr); ERFA._free(dcPtr);
  }
}
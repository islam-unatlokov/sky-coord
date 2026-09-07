import * as Astronomy from 'https://cdn.jsdelivr.net/npm/astronomy-engine@2.1.19/+esm';
import { IERS_DATA } from 'https://cdn.jsdelivr.net/gh/islam-unatlokov/sky-coord@main/iers.js';

console.log(IERS_DATA.updated_at);

const DEG2RAD = Math.PI / 180.0;
const RAD2DEG = 180.0 / Math.PI;
const ARCSEC2RAD = DEG2RAD / 3600.0;
const SPEED_OF_LIGHT_M_S = 299792458.0;
const SPEED_OF_LIGHT_AU_PER_DAY = 173.1446326847;
const OMEGA_EARTH = 7.292115146706979e-5; // rad/s (IERS nominal Earth rotation rate)

// WGS84 Ellipsoid Constants
const WGS84_A = 6378137.0; 
const WGS84_F = 1.0 / 298.257223563;
const WGS84_E2 = 2 * WGS84_F - WGS84_F * WGS84_F;

let iersData = null;

/**
 * Loads IERS Bulletin A table for UT1 - UTC and polar motion corrections.
 */
export async function loadIERSData() {
  if (!IERS_DATA) {
    throw new Error('[IERS Error] Could not find "IERS_DATA". Ensure iers.js is loaded in index.html.');
  }

  if (!IERS_DATA.data) {
    throw new Error('[IERS Error] Invalid format in "IERS_DATA". Missing "data" key.');
  }

  iersData = IERS_DATA.data;
  // return iersData;
}

/**
 * Returns interpolated IERS parameters and reports fallback status.
 */
export function getIERS(unixTime) {
  if (!iersData) {
    console.warn("[IERS Warning] Table not loaded. Using zero corrections (dUT1=0, xp=0, yp=0).");
    return { dUT1: 0.0, xp: 0.0, yp: 0.0, isExact: false };
  }

  const exactMJD = (unixTime / 86400) + 40587.0;
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
 * Core Coordinate Transformation Engine.
 * Converts topocentric horizontal coordinates to equatorial ICRS (RA, Dec).
 */
export function getEquatorialAstronomyEngine(lat_deg, lon_deg, elevation, zenith_deg, azimuth_deg, unix_time, { strict = true } = {}) {
  if (isNaN(azimuth_deg) || isNaN(zenith_deg) || isNaN(unix_time)) {
    throw new TypeError(`[Input Error] Received NaN: az=${azimuth_deg}, zen=${zenith_deg}, time=${unix_time}`);
  }

  const iers = getIERS(unix_time);
  if (strict && !iers.isExact) {
    throw new Error(`[Precision Error] Missing valid IERS parameters for timestamp ${unix_time}.`);
  }

  const { dUT1, xp, yp } = iers;

  // --- STEP 1: Topocentric Direction Vector (East, North, Up) ---
  const azRad = azimuth_deg * DEG2RAD;
  const zenRad = zenith_deg * DEG2RAD;
  const altRad = (90.0 - zenith_deg) * DEG2RAD;

  const xLocal = Math.cos(altRad) * Math.sin(azRad); // East
  const yLocal = Math.cos(altRad) * Math.cos(azRad); // North
  const zLocal = Math.sin(altRad);                  // Up / Zenith

  // --- STEP 2: Topocentric Frame to ITRS (WGS84 Normal) ---
  const latRad = lat_deg * DEG2RAD;
  const lonRad = lon_deg * DEG2RAD;

  const sinLat = Math.sin(latRad);
  const cosLat = Math.cos(latRad);
  const sinLon = Math.sin(lonRad);
  const cosLon = Math.cos(lonRad);

  // WGS84 Topocentric Basis Vectors in ITRS
  const E = [-sinLon, cosLon, 0.0];
  const N = [-sinLat * cosLon, -sinLat * sinLon, cosLat];
  const U = [ cosLat * cosLon,  cosLat * sinLon, sinLat];

  const rITRS = [
    xLocal * E[0] + yLocal * N[0] + zLocal * U[0],
    xLocal * E[1] + yLocal * N[1] + zLocal * U[1],
    xLocal * E[2] + yLocal * N[2] + zLocal * U[2]
  ];

  // --- STEP 3: Apply Correct Polar Motion W(xp, yp) Matrix -> TIRS ---
  const xpRad = xp * ARCSEC2RAD;
  const ypRad = yp * ARCSEC2RAD;

  // SOFA W Matrix: W13 = -xp, W23 = +yp, W31 = +xp, W32 = -yp
  const rTIRS = [
    rITRS[0] - xpRad * rITRS[2],
    rITRS[1] + ypRad * rITRS[2],
    xpRad * rITRS[0] - ypRad * rITRS[1] + rITRS[2]
  ];

  // --- STEP 4: Remove Diurnal Aberration ---
  const N_phi = WGS84_A / Math.sqrt(1.0 - WGS84_E2 * sinLat * sinLat);
  const xObs = (N_phi + elevation) * cosLat * cosLon;
  const yObs = (N_phi + elevation) * cosLat * sinLon;

  // Observer velocity in ITRF due to Earth rotation: v = w x r
  const vx_diurnal = -OMEGA_EARTH * yObs;
  const vy_diurnal =  OMEGA_EARTH * xObs;

  const beta_diurnal_x = vx_diurnal / SPEED_OF_LIGHT_M_S;
  const beta_diurnal_y = vy_diurnal / SPEED_OF_LIGHT_M_S;

  // Remove diurnal aberration from topocentric direction
  const rGeocentric_TIRS = [
    rTIRS[0] - beta_diurnal_x,
    rTIRS[1] - beta_diurnal_y,
    rTIRS[2]
  ];

  const normGeocentric = Math.hypot(rGeocentric_TIRS[0], rGeocentric_TIRS[1], rGeocentric_TIRS[2]);
  const vTIRS = [
    rGeocentric_TIRS[0] / normGeocentric,
    rGeocentric_TIRS[1] / normGeocentric,
    rGeocentric_TIRS[2] / normGeocentric
  ];

  // --- STEP 5: Rotate TIRS to True Equator of Date (EQD) via GAST at UT1 ---
  const timeUT1 = Astronomy.MakeTime(new Date((unix_time + dUT1) * 1000));
  const gastHours = Astronomy.SiderealTime(timeUT1);
  const gastRad = gastHours * 15.0 * DEG2RAD;

  const cosG = Math.cos(gastRad);
  const sinG = Math.sin(gastRad);

  const vecEQD = {
    x: vTIRS[0] * cosG - vTIRS[1] * sinG,
    y: vTIRS[0] * sinG + vTIRS[1] * cosG,
    z: vTIRS[2]
  };

  // --- STEP 6: Precession / Nutation (EQD -> EQJ / J2000) ---
  const timeUTC = Astronomy.MakeTime(new Date(unix_time * 1000));
  const rot = Astronomy.Rotation_EQD_EQJ(timeUTC);
  const vecJ2000 = Astronomy.RotateVector(rot, vecEQD);

  // --- STEP 7: Remove Annual Relativistic Aberration (SOFA eraAb Inverse) ---
  const earthState = Astronomy.BaryState(Astronomy.Body.Earth, timeUTC);
  const beta_x = earthState.vx / SPEED_OF_LIGHT_AU_PER_DAY;
  const beta_y = earthState.vy / SPEED_OF_LIGHT_AU_PER_DAY;
  const beta_z = earthState.vz / SPEED_OF_LIGHT_AU_PER_DAY;

  const px = vecJ2000.x;
  const py = vecJ2000.y;
  const pz = vecJ2000.z;

  const p_dot_v = px * beta_x + py * beta_y + pz * beta_z;

  // Exact first/second-order inverse aberration projection
  const p0_x = px - beta_x + p_dot_v * px;
  const p0_y = py - beta_y + p_dot_v * py;
  const p0_z = pz - beta_z + p_dot_v * pz;

  const normP0 = Math.hypot(p0_x, p0_y, p0_z);
  const vecAstrometric = {
    x: p0_x / normP0,
    y: p0_y / normP0,
    z: p0_z / normP0
  };

  // --- STEP 8: Convert to ICRS RA / Dec ---
  const sphere = Astronomy.SphereFromVector(vecAstrometric);
  let ra_deg = sphere.lon;
  if (ra_deg < 0) ra_deg += 360.0;
  const dec_deg = sphere.lat;

  return { ra_deg, dec_deg };
}

// Auto-initialize IERS table on module load
try {
  await loadIERSData();
} catch (err) {
  console.error("IERS Initialization Failed:", err);
}

// // Auto-initialize IERS data loading
// if (typeof document !== 'undefined') {
//   if (document.readyState === 'loading') {
//     document.addEventListener('DOMContentLoaded', () => loadIERSData());
//   } else {
//     loadIERSData();
//   }
// } else {
//   loadIERSData();
// }

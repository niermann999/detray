/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "detray/core/detector.hpp"

#include <iostream>

namespace detray {


/*enum mask_id : unsigned int {
    e_mask_types = 5,
    e_rectangle2 = 0,
    e_trapezoid2 = 1,
    e_annulus2 = 2,
    e_cylinder3 = 3,
    e_portal_cylinder3 = 3,  // no distinction from surface cylinder
    e_portal_ring2 = 4,
    e_single3 = std::numeric_limits<unsigned int>::max(),
    e_unknown = std::numeric_limits<unsigned int>::max(),
};*/

/** Creates a number of pixel modules for the a cylindrical barrel region.
 *
 * @tparam surface_t The surface type that contains the container indices
 * @tparam mask_t The geomterical boundaries of a surface (rectangles)
 *
 * @param m_half_x module half length in local x
 * @param m_half_y module half length in local y
 * @param m_tilt_phi phi tilt of the modules
 * @param layer_r radius at which the modules are positioned in the volume
 * @param radial_stagger module stagger in r
 * @param l_overlap the z overlap of modules next in z
 * @param binning phi, z bins of the surface grid
 *
 * @return a tuple that contains the surfaces (linking into the locally
 *         created container), the module transsforms and the surface masks.
 */
template <typename surface_t, typename mask_t>
inline auto create_modules(const scalar m_half_x = 8.4,
                           const scalar m_half_y = 36.,
                           const scalar m_tilt_phi = 0.14,
                           //const scalar m_tilt_phi = 0.145,
                           const scalar layer_r = 32.,
                           const scalar radial_stagger = 0.5,
                           const scalar l_overlap = 2.,
                           //const scalar radial_stagger = 2.,
                           //const scalar l_overlap = 5.,
                           const std::pair<int, int> binning = {16, 14}) {

    // Algebra type definitions from the plugins
    using point3 = __plugin::point3;
    using vector3 = __plugin::vector3;
    using transform3 = __plugin::transform3;

    /// mask index: type, range
    using mask_index = detray::darray<detray::dindex, 2>;

    // Prepare the return container
    using surface_container = detray::dvector<surface_t>;
    using trfs_container = detray::dvector<transform3>;
    using mask_container = detray::dvector<mask_t>;

    surface_container surfaces;
    trfs_container transforms;
    mask_container masks;

    // Create the module centers

    // surface grid bins
    int n_phi_bins = binning.first;
    int n_z_bins = binning.second;
    // module positions
    detray::dvector<point3> m_centers;
    m_centers.reserve(n_phi_bins * n_z_bins);

    // prep work
    scalar pi{static_cast<scalar>(M_PI)};
    scalar phi_step = scalar{2} * pi / (n_phi_bins);
    scalar min_phi = -pi + scalar{0.5} * phi_step;
    scalar z_start =
        scalar{-0.5} * (n_z_bins - 1) * (scalar{2} * m_half_y - l_overlap);
    scalar z_step = scalar{2} * std::abs(z_start) / (n_z_bins - 1);

    // loop over the z bins
    for (size_t z_bin = 0; z_bin < size_t(n_z_bins); ++z_bin) {
        // prepare z and r
        scalar m_z = z_start + z_bin * z_step;
        scalar m_r = (z_bin % 2) != 0u ? layer_r - scalar{0.5} * radial_stagger
                                       : layer_r + scalar{0.5} * radial_stagger;
        for (size_t phiBin = 0; phiBin < size_t(n_phi_bins); ++phiBin) {
            // calculate the current phi value
            scalar m_phi = min_phi + phiBin * phi_step;
            m_centers.push_back(
                point3{m_r * std::cos(m_phi), m_r * std::sin(m_phi), m_z});
        }
    }

    // Create geometry data

    // First value is two for rectangle type, then index into local container
    mask_index m_id = {0, 0};

    for (auto& m_center : m_centers) {

        // Surfaces with the linking into the local containers
        m_id = {0, masks.size()};
        surfaces.emplace_back(transforms.size(), m_id, detray::dindex_invalid,
                              detray::dindex_invalid);

        // The rectangle bounds for this module
        masks.emplace_back(m_half_x, m_half_y);
        masks.back().links() = {dindex_invalid, dindex_invalid};

        // Build the transform
        // The local phi
        scalar m_phi = algebra::getter::phi(m_center);
        // Local z axis is the normal vector
        vector3 m_local_z{std::cos(m_phi + m_tilt_phi),
                          std::sin(m_phi + m_tilt_phi), 0.};
        // Local x axis the normal to local y,z
        vector3 m_local_x{-std::sin(m_phi + m_tilt_phi),
                          std::cos(m_phi + m_tilt_phi), 0.};

        // Create the module transform
        transforms.emplace_back(m_center, m_local_z, m_local_x);
    }

    return std::make_tuple<surface_container, trfs_container, mask_container>(
        std::move(surfaces), std::move(transforms), std::move(masks));
}

/** Helper method for positioning
  *
  * @param z is the z position of the ring
  * @param radius is the ring radius
  * @param phi_stagger is the radial staggering along phi
  * @param phi_sub_stagger is the overlap of the modules
  * @param n_phi_bins is the number of bins in phi
  *
  * @return a vector of the module positions in a ring
  */
inline auto module_positions_ring(scalar z,
                                   scalar radius,
                                   scalar phi_stagger,
                                   scalar phi_sub_stagger,
                                   int n_phi_bins) {
    // create and fill the positions
    std::vector<vector3> r_positions;
    r_positions.reserve(n_phi_bins);

    // prep work
    scalar pi{static_cast<scalar>(M_PI)};
    scalar phi_step = scalar{2} * pi / (n_phi_bins);
    double min_phi = -pi + 0.5 * phi_step;

    for (size_t iphi = 0; iphi < size_t(n_phi_bins); ++iphi) {
        // if we have a phi sub stagger presents
        double rzs = 0.;
        // phi stagger affects 0 vs 1, 2 vs 3 ... etc
        // -> only works if it is a %4
        // phi sub stagger affects 2 vs 4, 1 vs 3 etc.
        if (phi_sub_stagger != 0. && !(n_phi_bins % 4)) {
            // switch sides
            if (!(iphi % 4)) {
                rzs = phi_sub_stagger;
            }
            else if (!((iphi + 1) % 4)) {
                rzs = -phi_sub_stagger;
            }
        }
        // the module phi
        double phi = min_phi + iphi * phi_step;
        // main z position depending on phi bin
        double rz = iphi % 2 ? z - 0.5 * phi_stagger : z + 0.5 * phi_stagger;
        r_positions.push_back(
            vector3{radius * std::cos(phi), radius * std::sin(phi), rz + rzs});
    }
    return r_positions;
}

/** Helper method for positioning rings in a disc
*
* @param z is the nominal z posiiton of the dis
* @param ring_stagger is the staggering of the different rings
* @param phi_stagger is the staggering on a ring in phi : it is even/odd
* @param phi_sub_stagger is the sub staggering on a ring in phi : it affects
* 0/4/8 and 3/6
* @param inner_r is the inner Radius for the disc
* @param outer_r is the outer Radius for the disc
* @param disc_binning is the binning setup in r (size of the vector), phi
* @param m_half_y is a pair of phibins and module length
* @param m_half_x_min_y The half lenght in X (at Y min) of the module
* @param m_half_x_max_y The half lenght in X (at Y max) of the module
* @param m_half_y The half lenght in Y of the module
* @param m_tilt The tilt out of the plane for discs
*
* @return vector of module positions of a ring
*/
template <typename surface_t, typename mask_t>
inline auto create_endcap_modules(scalar z,
                                  scalar ring_stagger = 0.0, 
                                  std::vector<scalar> phi_stagger = {4.0, 4.0},
                                  std::vector<scalar> phi_sub_stagger = {0.5, 0.},
                                  scalar inner_r = 27.,
                                  scalar outer_r = 180.,
                                  const std::vector<size_t>& disc_binning = {40, 68},
                                  const std::vector<scalar>& m_half_y = {36., 36.},
                                  std::vector<scalar> m_half_x_min_y = {8.4, 8.4}, 
                                  std::vector<scalar> m_half_x_max_y = {12.4, 12.4}, 
                                  std::vector<scalar> m_tilt = {0., 0.},
                                  int side = 1) {
    // Algebra type definitions from the plugins
    using point3 = __plugin::point3;
    using vector3 = __plugin::vector3;
    using transform3 = __plugin::transform3;

    /// mask index: type, range
    using mask_index = detray::darray<detray::dindex, 2>;

    // Prepare the return container
    using surface_container = detray::dvector<surface_t>;
    using trfs_container = detray::dvector<transform3>;
    using mask_container = detray::dvector<mask_t>;

    surface_container surfaces;
    trfs_container transforms;
    mask_container masks;


    // calculate the radii of the rings
    std::vector<scalar> radii;
    // calculate the radial borders
    std::vector<scalar> radial_boarders;
    // the radial span of the disc
    scalar delta_r = outer_r - inner_r;

    // Only one ring
    if (disc_binning.size() == 1) {
        radii.push_back(scalar{0.5} * (inner_r + outer_r));
        radial_boarders = {inner_r, outer_r};
    }
    else {
        // sum up the total length of the modules along r
        scalar tot_length = 0;
        for (auto& mhlength : m_half_y) {
            tot_length += scalar{2} * mhlength;
        }
        // now calculate the overlap (equal pay)
        scalar r_overlap = (tot_length - delta_r) / (m_half_y.size() - 1);
        // and now fill the radii and gaps
        scalar last_r = inner_r;
        scalar last_hl = 0.;
        scalar last_ol = 0.;
        // remember the radial boarders
        radial_boarders.push_back(inner_r);
        for (auto& mhlength : m_half_y) {
            // calculate the radius
            radii.push_back(last_r + last_hl - last_ol + mhlength);
            last_r = radii.back();
            last_ol = r_overlap;
            last_hl = mhlength;
            // and register the radial boarder
            radial_boarders.push_back(last_r + scalar{2} * last_hl - scalar{0.5} * last_ol);
        }
    }

    // now build the modules in every ring

    // Mask type is not yet known, but its the first in its container
    mask_index m_id = {1, 0};

    for (size_t ir = 0; ir < radii.size(); ++ir) {
        // generate the z value
        // convention inner ring is closer to origin : makes sense
        double rz = radii.size() == 1
                        ? z
                        : (ir % 2 ? z + scalar{0.5} * ring_stagger : z - scalar{0.5} * ring_stagger);
        // fill the ring module positions
        double ps_stagger = phi_sub_stagger.size() ? phi_sub_stagger[ir] : 0.;
        /*_positions.push_back(module_positions_ring(rz, radii[ir], 
                                                    phi_stagger[ir],
                                                    ps_stagger,
                                                    disc_binning[ir]));*/

        std::vector<point3> r_postitions = module_positions_ring(rz, radii[ir], 
                                                    phi_stagger[ir],
                                                    ps_stagger,
                                                    disc_binning[ir]);

        // Build the geometrical objects
        for (const auto &m_position : r_postitions) {
            // rectangle
            /*if (m_half_x_min_y[ir] == m_half_x_max_y[ir]) {
                m_id = {0, masks.size()};
                masks.emplace_back(m_half_x_min_y[ir], m_half_y[ir]);
            // trapezoid
            } else {*/
            m_id = {1, masks.size()};
            masks.emplace_back(m_half_x_min_y[ir], m_half_x_max_y[ir], m_half_y[ir]);
            //}
            // The links will be updated to the volume later
            masks.back().links() = {dindex_invalid, dindex_invalid};

            // Surfaces with the linking into the local containers
            surfaces.emplace_back(transforms.size(), m_id, detray::dindex_invalid, detray::dindex_invalid);

            // the module transform from the position
            double m_phi = algebra::getter::phi(m_position);
            // the center position of the modules
            point3 m_center{static_cast<scalar>(side) * m_position};
            // the rotation matrix of the module
            vector3 m_local_y{std::cos(m_phi), std::sin(m_phi), 0.};
            // take different axis to have the same readout direction
            vector3 m_local_z{0., 0., side * 1.};
            vector3 m_local_x = algebra::vector::cross(m_local_y, m_local_z);

            // Create the module transform
            transforms.emplace_back(m_center, m_local_z, m_local_x);
        }
    }

    return std::make_tuple<surface_container, trfs_container, mask_container>(
        std::move(surfaces), std::move(transforms), std::move(masks));
}




/** Generate surfaces of a disc
  *
  * @param m_half_x_min_y The half lenght in X (at Y min) of the module
  * @param m_half_x_max_y The half lenght in X (at Y max) of the module
  * @param m_half_y The half lenght in Y of the module
  * @param moduleTilt The tilt out of the plane for discs
  * @param ringRadius The central radius of the ring
  * @param ringZ The z position of the ring
  * @param zStagger The z offset of phi moudles
  * @param nPhi The number of phi modules
  *
  * @return A vector of surfaces
  */

    /*ActsScalar moduleHalfXminY = 8.4, ActsScalar moudleHalfXmaxY = 12.4,
    ActsScalar moduleHalfY = 32., ActsScalar moduleTilt = 0.,
    ActsScalar ringRadius = 40., ActsScalar zStagger = 2, int nPhi = 40*/
/*template<typename rectangle_t, template trapezoid_t>
inline auto create_endcap_modules(
    scalar m_half_x_min_y = 8.4, scalar m_half_x_max_y = 12.4,
    scalar m_half_y = 32., scalar m_tiltt = 0., scalar ring_r = 40.,
    scalar ring_z = 0., scalar z_stagger = 2,, int n_phi = 40) {

    // Algebra type definitions from the plugins
    using point3 = __plugin::point3;
    using vector3 = __plugin::vector3;
    using transform3 = __plugin::transform3;

    /// mask index: type, range
    using mask_index = detray::darray<detray::dindex, 2>;

    // Prepare the return container
    using surface_container = detray::dvector<surface_t>;
    using trfs_container = detray::dvector<transform3>;
    using mask_container = detray::dvector<mask_t>;

    surface_container surfaces;
    trfs_container transforms;
    mask_container masks;

    // Prepare the return vector
    //std::vector<std::shared_ptr<Surface>> layerSurfaces;

    // The rectangle/trapezoid bounds for all modules
    //std::shared_ptr<PlanarBounds> mBounds = nullptr;


    // Mask type is not yet known, but its the first in its container
    mask_index m_id = {dindex_invalid, 0};

    // prep work
    scalar pi{static_cast<scalar>(M_PI)};
    scalar phi_step = scalar{2} * pi / (n_phi);

    for (int im = 0; im < nPhi; ++im) {
        // Get the moduleTransform
        scalar phi = -pi + im * phi_step;
        auto mModuleTransform = Transform3(
            Translation3(ringRadius * std::cos(phi), ringRadius * std::sin(phi),
                        ringZ + (im % 2) * zStagger) *
            AngleAxis3(phi - 0.5 * M_PI, Vector3::UnitZ()) *
            AngleAxis3(moduleTilt, Vector3::UnitY()));

        // Create the detector element
        //auto detSurface =
        //    Surface::makeShared<PlaneSurface>(mModuleTransform, mBounds);
        //layerSurfaces.push_back(detSurface);

        // rectangle
        if (m_half_x_min_y == m_half_x_max_y) {
            m_id = {2, masks.size()};
            masks.emplace_back(m_half_x_min_y, m_half_y);
        // trapezoid
        } else {
            m_id = {3, masks.size()};
            masks.emplace_back(m_half_x_min_y, m_half_x_max_y, m_half_y);
        }
        // The rectangle bounds for this module
        masks.back().links() = {dindex_invalid, dindex_invalid};

        // Surfaces with the linking into the local containers
        surfaces.emplace_back(transforms.size(), m_id, detray::dindex_invalid,
                              detray::dindex_invalid);


     // low loop over the phi positions and build the stuff
    for (auto& ringModulePosition : discModulePositions) {
        // the module transform from the position
        double m_phi = algebra::getter::phi(ringModulePosition);
        // the center position of the modules
        point3 m_center{side * ringModulePosition};
        // the rotation matrix of the module
        vector3 m_local_y{std::cos(m_phi), std::sin(m_phi), 0.};
        // take different axis to have the same readout direction
        vector3 m_local_z{0., 0., side * 1.};
        vector3 m_local_x = algebra::vector::cross(m_local_y, m_local_z);

        // Create the module transform
        transforms.emplace_back(m_center, m_local_z, m_local_x);
    }
    return layerSurfaces;
}*/

/// @return an endcap volume
/*std::shared_ptr<DetectorVolume> createEndcapVolume(
    scalar volumeMinR, scalar volumeMaxR, scalar volumeMinZ,
    scalar volumeMaxZ, const std::string& volumeName = "SingleLayerVolume",
    scalar moduleHalfXminY = 8.4, scalar moudleHalfXmaxY = 12.4,
    scalar moduleHalfY = 32., scalar moduleTilt = 0.,
    scalar ringRadius = 40., scalar zStagger = 2, int nPhi = 40) {
  // Place the ring into the middle
  scalar ringZ = 0.5 * (volumeMinZ + volumeMaxZ);

  auto volumeSurfaces =
      surfacesRing(moduleHalfXminY, moudleHalfXmaxY, moduleHalfY, moduleTilt,
                   ringRadius, ringZ, zStagger, nPhi);

  // Create the volume bounds
  auto volumeBounds = std::make_unique<CylinderVolumeBounds>(
      volumeMinR, volumeMaxR, 0.5 * (volumeMaxZ - volumeMinZ));

  SurfaceLinks volumeSurfaceLinks = AllSurfaces{};

  std::vector<SurfaceLinks> portalSurfaceLinks = {AllSurfaces{}, AllSurfaces{},
                                                  AllSurfaces{}};

  if (volumeMinR > 0.) {
    portalSurfaceLinks.push_back(AllSurfaces{});
  }

  auto volumeTransform = Transform3::Identity();
  volumeTransform.pretranslate(Vector3(0., 0., ringZ));

  return DetectorVolume::makeShared(
      volumeTransform, std::move(volumeBounds), std::move(volumeSurfaces),
      std::move(volumeSurfaceLinks), std::move(portalSurfaceLinks), volumeName);
}*/




/// Helper method to create a central detector
///
/// Keep track of the central half length with
/// @param detectorRmin inner radius of detector
/// @param detectorRmax outer radius of detector
/// @param zToCentral is the distance to central to
/// @param side is the side of the endcap detector
/// @param detectorName is the detector name prescript
///
/// @return a central detector volume
/*std::shared_ptr<DetectorVolume> createEndcapDetector(
    scalar detectorRmin = 0., scalar detectorRmax = 80.,
    scalar zToCentral = 500., int side = 1) {
  scalar layerThickness = 5.;
  scalar gapThickness = 50.;

  std::string firstLayerName = (side > 0) ? "Layer0" : "Layer1";
  std::string gapName = "Gap";
  std::string secondLayerName = (side > 0) ? "Layer1" : "Layer0";
  std::string sideTag = (side > 0) ? "Pos" : "Neg";

  // Place the first layer
  scalar oneZ = side * (zToCentral);
  scalar twoZ = side * (zToCentral + layerThickness);
  auto firstLayer = createEndcapVolume(
      detectorRmin, detectorRmax, std::min(oneZ, twoZ), std::max(oneZ, twoZ),
      detectorName + firstLayerName + sideTag);

  // Adapt for the gap & build
  oneZ = side * (zToCentral + layerThickness);
  twoZ = side * (zToCentral + layerThickness + gapThickness);
  Transform3 gapTransform = Transform3::Identity();
  gapTransform.pretranslate(Vector3(0., 0., 0.5 * (oneZ + twoZ)));
  auto gapBounds = std::make_unique<CylinderVolumeBounds>(
      detectorRmin, detectorRmax, std::abs(0.5 * (oneZ - twoZ)));
  auto gap = DetectorVolume::makeShared(gapTransform, std::move(gapBounds),
                                        detectorName + gapName + sideTag);

  // Adapt for the second layer
  oneZ = side * (zToCentral + layerThickness + gapThickness);
  twoZ = side * (zToCentral + 2 * layerThickness + gapThickness);
  auto secondLayer = createEndcapVolume(
      detectorRmin, detectorRmax, std::min(oneZ, twoZ), std::max(oneZ, twoZ),
      detectorName + secondLayerName + sideTag);

  std::vector<std::shared_ptr<DetectorVolume>> endcapVolumes;
  if (side > 0) {
    endcapVolumes = {firstLayer, gap, secondLayer};
  } else {
    endcapVolumes = {secondLayer, gap, firstLayer};
  }
  // Container in Z
  return CylindricalContainerHelper::containerInZ(
      std::move(endcapVolumes),
      detectorName + std::string("TwoLayers") + sideTag);
}*/

/** Builds a simple detray geometry of the innermost tml layers. It contains:
 *
 * - a beampipe (r = 27mm, half_z = 500mm)
 * - a first layer (r_min = 27mm, r_max = 38mm, half_z = 500mm) with 224
 *   rectangular (half_x = 8.4mm, half_y = 36mm) modules at r = 32mm
 * - an empty layer (r_min = 38mm, r_max = 64mm, half_z = 500mm)
 * - a second layer (r_min = 64mm, r_max = 80mm, half_z = 500mm) with 448
 *   rectangular (half_x = 8.4mm, half_y = 36mm) modules at r = 72mm.
 *
 * @returns a tuple containing the geometry objects collections: [volumes,
 *          surfaces, transforms, disc masks (neg/pos portals), cylinder masks
 *          (inner/outer portals), rectangle masks (modules)]
 */
auto create_toy_geometry(vecmem::memory_resource& resource) {

    using transform3 = __plugin::transform3;

    // Volume type
    using volume_type =
        detray::volume<toy_object_registry, dindex_range, detray::darray>;
    /// volume index: volume the surface belongs to
    using volume_index = detray::dindex;
    /// transform link: transform entry belonging to surface
    using transform_link = detray::dindex;
    /// mask index: type, range
    using mask_index = detray::darray<detray::dindex, 2>;
    /// volume links: next volume, next (local) object finder
    using edge_links = detray::darray<detray::dindex, 2>;
    /// source link
    using source_link = detray::dindex;

    // We have rectangle and trapezoid surfaces for modules as well as discs 
    // (ring) and cylinder surfaces for portals
    using rectangle = detray::rectangle2<detray::planar_intersector,
                                         __plugin::cartesian2, edge_links, 0>;

    using trapezoid = detray::trapezoid2<detray::planar_intersector,
                                         __plugin::cartesian2, edge_links, 1>;

    using annulus = detray::annulus2<detray::planar_intersector,
                                         __plugin::cartesian2, edge_links, 2>;

    using cylinder = detray::cylinder3<false, detray::cylinder_intersector,
                                       __plugin::cylindrical2, edge_links, 3>;

    using disc = detray::ring2<detray::planar_intersector, __plugin::cartesian2,
                               edge_links, 4>;

    // The surface type, both for volume portals and contained detector
    // surfaces
    using surface = detray::surface_base<transform_link, mask_index,
                                         volume_index, source_link, edge_links>;

    // The geometry data containers
    using volume_container = detray::dvector<volume_type>;
    using surface_container = detray::dvector<surface>;
    using transf_container = detray::dvector<transform3>;
    using rectangle_container = detray::dvector<rectangle>;
    using trapezoid_container = detray::dvector<trapezoid>;
    using annulus_container = detray::dvector<annulus>;
    using cylinder_container = detray::dvector<cylinder>;
    using disc_container = detray::dvector<disc>;

    /** Volumes */
    volume_container volumes = {};
    /** Surfaces, including portals */
    surface_container surfaces = {};
    /** Surface transforms */
    transf_container transforms = {};
    /** Rectangle masks for detector barrel surfaces */
    rectangle_container rectangles = {};
    /** Trapezoid masks for detector endcap surfaces */
    trapezoid_container trapezoids = {};
    /** Empty*/
    annulus_container annuli = {};
    /** Cylinder masks for inner/outer volume boundaries (portals) */
    cylinder_container cylinders = {};
    /** Disc masks for neg/pos volume boundaries (portals) */
    disc_container discs = {};
    // mask index for surfaces
    mask_index m_id = {};
    /** source link */
    const dindex inv_sf_finder = dindex_invalid;
    /** Leaving world */
    const dindex leaving_world = dindex_invalid;

    //
    // general
    //
    const scalar detectorRmin = 0.;
    const scalar detectorRmax = 180.;
    const scalar beampipe_r = 27.;

    //
    // barrel
    //
    const scalar brl_half_z = 500.;
    const std::vector<scalar> brl_positions = {19., 32., 72., 116., 172.};
    const std::vector<std::pair<scalar, scalar>> brl_vol_sizes = {{0., 27.}, {27., 38.}, {64., 80.}, {108., 124.}, {164., 180.}};
    const scalar brl_radial_stagger = 0.5;
    const scalar brl_l_overlap = 2.;
    //const scalar brl_radial_stagger = 2.;
    //const scalar brl_l_overlap = 5.;
    const std::vector<std::pair<int, int>> brl_binning = {{0., 0.}, {16, 14}, {32, 14}, {52, 14}, {78, 14}};
    // module parameters
    const scalar brl_half_x = 8.4;
    const scalar brl_half_y = 36.;
    const scalar brl_tilt_phi = 0.14;
    //const scalar brl_tilt_phi = 0.145;

    //
    // endcaps
    //
    const std::vector<scalar> edc_positions = {600., 700., 820., 960., 1100., 1300., 1500.};
    const std::vector<std::pair<scalar, scalar>> edc_vol_sizes = {{595., 605.}, {695., 705.}, {815., 825.}, {955., 965.}, {1095., 1105.}, {1295., 1305.}, {1495., 1505.}};
    const scalar edc_ring_stagger = 0.0;
    // Parameters for both rings of modules
    const std::vector<scalar> edc_phi_stagger = {4.0, 4.0};
    const std::vector<scalar> edc_phi_sub_stagger = {0.5, 0.};
    const scalar edc_inner_r = 27.;
    const scalar edc_outer_r = 180.;
    const std::vector<size_t>& edc_disc_binning = {40, 68};
    // module params
    const std::vector<scalar>& edc_half_y = {36., 36.};
    const std::vector<scalar> edc_half_x_min_y = {8.4, 8.4};
    const std::vector<scalar> edc_half_x_max_y = {12.4, 12.4};
    const std::vector<scalar> edc_tilt = {0., 0.};

    struct toy_grid {};

    /** Add a single barrel layer volume to an existing collection.
     *
     * @param min_r minimal radius of volume
     * @param max_r maximal radius of volume
     * @param half_z half length in z of volume
     *
     * @return a detray cylinder volume
     */
    auto add_cylinder_volume = [&](const scalar min_r, const scalar max_r,
                                   const scalar lower_z, const scalar upper_z) -> volume_type& {
        // The volume bounds
        detray::darray<scalar, 6> bounds = {min_r,  max_r,
                                            std::min(lower_z, upper_z),
                                            std::max(lower_z, upper_z),
                                            -M_PI, M_PI};
        // Add the new volume to the collection
        auto& vol = volumes.emplace_back(bounds);
        vol.set_index(volumes.size() - 1);

        return vol;
    };

    /** Function that adds a cylinder portal.
     *
     * @param vol volume the portal should be added to
     * @param r raius of the cylinder
     * @param half_z z half length of the cylinder
     */
    auto add_cylinder_portal = [&](volume_type& vol, const scalar r,
                                   const scalar lower_z, const scalar upper_z) {
        m_id = {3, cylinders.size()};
        surfaces.emplace_back(transforms.size(), m_id, vol.index(),
                              inv_sf_finder);
        cylinders.emplace_back(r, std::min(lower_z, upper_z), 
                               std::max(lower_z, upper_z));
        cylinders.back().links() = {dindex_invalid, dindex_invalid};
        // Position cylinder
        //vector3 translation{0., 0., scalar{0.5} * (upper_z + lower_z)};
        vector3 translation{0., 0., 0.};
        transforms.emplace_back(translation);

        vol.update_range({surfaces.size() - 1, surfaces.size()});
    };

    /** Function that adds a disc portal.
     *
     * @param vol volume the portal should be added to
     * @param min_r lower radius of disc
     * @param max_r upper radius of disc
     * @param half_z z half length of the detector volume
     */
    auto add_disc_portal = [&](volume_type& vol, const scalar min_r,
                               const scalar max_r, const scalar z) {
        m_id = {4, discs.size()};
        surfaces.emplace_back(transforms.size(), m_id, vol.index(),
                              inv_sf_finder);
        discs.emplace_back(min_r, max_r);
        discs.back().links() = {dindex_invalid, dindex_invalid};
        // Position disc
        vector3 translation{0., 0., z};
        transforms.emplace_back(translation);

        vol.update_range({surfaces.size() - 1, surfaces.size()});
    };

    /** Function that updates surface and volume links when added to the global
     *  containers. The module surface edge point back into the volume.
     *
     * @param vol the volume the module surfaces belong to
     * @param modules the module surface container
     * @param trfs_offset offset into the global transform container
     * @param masks_offset offset into the global mask container
     */
    auto update_links = [&](volume_type& vol, surface_container& modules,
                            const dindex trfs_offset,
                            const dindex masks_offset) {
        for (auto& sf : modules) {
            sf.transform() += trfs_offset;
            std::get<1>(sf.mask()) += masks_offset;
            sf.volume() = vol.index();
            sf.set_edge({vol.index(), inv_sf_finder});
        }

        vol.update_range({surfaces.size(), surfaces.size() + modules.size()});
    };

    //
    // beampipe
    //

    // build the beampipe volume
    volume_type& beampipe =
        add_cylinder_volume(brl_vol_sizes[0].first, brl_vol_sizes[0].second, -edc_vol_sizes[2].second, edc_vol_sizes[2].second);

    // This is the beampipe surface
    m_id = {3, cylinders.size()};
    surfaces.emplace_back(transforms.size(), m_id, beampipe.index(),
                              inv_sf_finder);
    surfaces.back().set_edge({beampipe.index(), inv_sf_finder});
    cylinders.emplace_back(brl_positions[0], -edc_vol_sizes[2].second, edc_vol_sizes[2].second);
    cylinders.back().links() = {beampipe.index(), inv_sf_finder};
    transforms.emplace_back();  // identity
    beampipe.update_range({surfaces.size() - 1, surfaces.size()});

    // negative and positive, outer portal surfaces
    // cylinder portals for all volumes
    //negative endcap
    add_cylinder_portal(beampipe, edc_inner_r, -edc_vol_sizes[2].second, -edc_vol_sizes[2].first);
    add_cylinder_portal(beampipe, edc_inner_r, -edc_vol_sizes[2].first, -edc_vol_sizes[1].second);
    add_cylinder_portal(beampipe, edc_inner_r, -edc_vol_sizes[1].second, -edc_vol_sizes[1].first);
    add_cylinder_portal(beampipe, edc_inner_r, -edc_vol_sizes[1].first, -edc_vol_sizes[0].second);
    add_cylinder_portal(beampipe, edc_inner_r, -edc_vol_sizes[0].second, -edc_vol_sizes[0].first);
    add_cylinder_portal(beampipe, edc_inner_r, -edc_vol_sizes[0].first, -brl_half_z);
    // barrel
    add_cylinder_portal(beampipe, edc_inner_r, -brl_half_z, brl_half_z);
    // positive endcap
    add_cylinder_portal(beampipe, edc_inner_r, brl_half_z, edc_vol_sizes[0].first);
    add_cylinder_portal(beampipe, edc_inner_r, edc_vol_sizes[0].first, edc_vol_sizes[0].second);
    add_cylinder_portal(beampipe, edc_inner_r, edc_vol_sizes[0].second, edc_vol_sizes[1].first);
    add_cylinder_portal(beampipe, edc_inner_r, edc_vol_sizes[1].first, edc_vol_sizes[1].second);
    add_cylinder_portal(beampipe, edc_inner_r, edc_vol_sizes[1].second, edc_vol_sizes[2].first);
    add_cylinder_portal(beampipe, edc_inner_r, edc_vol_sizes[2].first, edc_vol_sizes[2].second);
    // discs
    add_disc_portal(beampipe, brl_vol_sizes[0].first, brl_vol_sizes[0].second, -edc_vol_sizes[2].second);
    add_disc_portal(beampipe, brl_vol_sizes[0].first, brl_vol_sizes[0].second, edc_vol_sizes[2].second);

    // Set disc portal edges
    dindex beampipe_idx = beampipe.index();
    dindex beampipe_pt_index = beampipe.range()[0] + 1;
    // negative endcap
    surfaces[beampipe_pt_index].set_edge({beampipe_idx + 1, inv_sf_finder});
    surfaces[++beampipe_pt_index].set_edge({beampipe_idx + 2, inv_sf_finder});
    surfaces[++beampipe_pt_index].set_edge({beampipe_idx + 3, inv_sf_finder});
    surfaces[++beampipe_pt_index].set_edge({beampipe_idx + 4, inv_sf_finder});
    surfaces[++beampipe_pt_index].set_edge({beampipe_idx + 5, inv_sf_finder});
    surfaces[++beampipe_pt_index].set_edge({beampipe_idx + 6, inv_sf_finder});
    // barrel
    surfaces[++beampipe_pt_index].set_edge({beampipe_idx + 7, inv_sf_finder});
    // positive endcap
    surfaces[++beampipe_pt_index].set_edge({beampipe_idx + 14, inv_sf_finder});
    surfaces[++beampipe_pt_index].set_edge({beampipe_idx + 15, inv_sf_finder});
    surfaces[++beampipe_pt_index].set_edge({beampipe_idx + 16, inv_sf_finder});
    surfaces[++beampipe_pt_index].set_edge({beampipe_idx + 17, inv_sf_finder});
    surfaces[++beampipe_pt_index].set_edge({beampipe_idx + 18, inv_sf_finder});
    surfaces[++beampipe_pt_index].set_edge({beampipe_idx + 19, inv_sf_finder});
    //discs
    surfaces[++beampipe_pt_index].set_edge({leaving_world, inv_sf_finder});
    surfaces[++beampipe_pt_index].set_edge({leaving_world, inv_sf_finder});

    //
    // negative endcap
    //

    int side = -1;

    // build the outermost volume
    scalar lay_neg_z = side * edc_vol_sizes[2].second;
    scalar lay_pos_z = side * edc_vol_sizes[2].first;

    volume_type& neg_edc_3 =
        add_cylinder_volume(edc_inner_r, edc_outer_r, lay_neg_z, lay_pos_z);

    // create disc module surfaces
    auto [neg_edc_3_modules, neg_edc_3_trfs, neg_edc_3_masks] = 
                    create_endcap_modules<surface, trapezoid>(
                                                       side * edc_positions[2],
                                                       edc_ring_stagger,
                                                       edc_phi_stagger,
                                                       edc_phi_sub_stagger,
                                                       edc_inner_r, 
                                                       edc_outer_r,
                                                       edc_disc_binning,
                                                       edc_half_y,
                                                       edc_half_x_min_y,
                                                       edc_half_x_max_y,
                                                       edc_tilt);

    // update linking in volumes and surfaces
    update_links(neg_edc_3, neg_edc_3_modules, transforms.size(), trapezoids.size());
    // Append to collections
    surfaces.insert(surfaces.end(), neg_edc_3_modules.begin(), neg_edc_3_modules.end());
    transforms.insert(transforms.end(), neg_edc_3_trfs.begin(), neg_edc_3_trfs.end());
    trapezoids.insert(trapezoids.end(), neg_edc_3_masks.begin(), neg_edc_3_masks.end());

    // negative and positive, inner and outer portal surface
    add_cylinder_portal(neg_edc_3, edc_inner_r, lay_neg_z, lay_pos_z);
    add_cylinder_portal(neg_edc_3, edc_outer_r, lay_neg_z, lay_pos_z);
    add_disc_portal(neg_edc_3, edc_inner_r, edc_outer_r, lay_neg_z);
    add_disc_portal(neg_edc_3, edc_inner_r, edc_outer_r, lay_pos_z);

    // Index of first portal
    dindex next_vol_idx = neg_edc_3.index() + 1;
    dindex layer_pt_index = neg_edc_3.range()[1] - 4;
    surfaces[layer_pt_index].set_edge({beampipe_idx, inv_sf_finder});
    surfaces[++layer_pt_index].set_edge({leaving_world, inv_sf_finder});
    surfaces[++layer_pt_index].set_edge({leaving_world, inv_sf_finder});
    surfaces[++layer_pt_index].set_edge({next_vol_idx, inv_sf_finder});

    //
    // gap layer
    //

    // build gap volume
    scalar gap_neg_z = side * edc_vol_sizes[2].first;
    scalar gap_pos_z = side * edc_vol_sizes[1].second;

    volume_type& neg_edc_gap_3 =
        add_cylinder_volume(edc_inner_r, edc_outer_r, gap_neg_z, gap_pos_z);

    // negative and positive, inner and outer portal surface
    add_cylinder_portal(neg_edc_gap_3, edc_inner_r,gap_neg_z, gap_pos_z);
    add_cylinder_portal(neg_edc_gap_3, edc_outer_r, gap_neg_z, gap_pos_z);
    add_disc_portal(neg_edc_gap_3, edc_inner_r, edc_outer_r, gap_neg_z);
    add_disc_portal(neg_edc_gap_3, edc_inner_r, edc_outer_r, gap_pos_z);

    // Connect portals (first layer, second layer)
    dindex prev_vol_idx = neg_edc_gap_3.index() - 1;
    next_vol_idx = neg_edc_gap_3.index() + 1;
    dindex gap_pt_index = neg_edc_gap_3.range()[1] - 4;
    // cylinder
    surfaces[gap_pt_index].set_edge({beampipe_idx, inv_sf_finder});
    surfaces[++gap_pt_index].set_edge({leaving_world, inv_sf_finder});
    // discs
    surfaces[++gap_pt_index].set_edge({prev_vol_idx, inv_sf_finder});
    surfaces[++gap_pt_index].set_edge({next_vol_idx, inv_sf_finder});

    //
    // second layer
    //
    lay_neg_z = side * edc_vol_sizes[1].second;
    lay_pos_z = side * edc_vol_sizes[1].first;

    volume_type& neg_edc_2 =
        add_cylinder_volume(edc_inner_r, edc_outer_r, lay_neg_z, lay_pos_z);

    // create disc module surfaces
    auto [neg_edc_2_modules, neg_edc_2_trfs, neg_edc_2_masks] = 
                    create_endcap_modules<surface, trapezoid>(
                                                       side * edc_positions[1],
                                                       edc_ring_stagger,
                                                       edc_phi_stagger,
                                                       edc_phi_sub_stagger,
                                                       edc_inner_r, 
                                                       edc_outer_r,
                                                       edc_disc_binning,
                                                       edc_half_y,
                                                       edc_half_x_min_y,
                                                       edc_half_x_max_y,
                                                       edc_tilt);

    // update linking in volumes and surfaces
    update_links(neg_edc_2, neg_edc_2_modules, transforms.size(), trapezoids.size());
    // Append to collections
    surfaces.insert(surfaces.end(), neg_edc_2_modules.begin(), neg_edc_2_modules.end());
    transforms.insert(transforms.end(), neg_edc_2_trfs.begin(), neg_edc_2_trfs.end());
    trapezoids.insert(trapezoids.end(), neg_edc_2_masks.begin(), neg_edc_2_masks.end());

    // negative and positive, inner and outer portal surface
    add_cylinder_portal(neg_edc_2, edc_inner_r, lay_neg_z, lay_pos_z);
    add_cylinder_portal(neg_edc_2, edc_outer_r, lay_neg_z, lay_pos_z);
    add_disc_portal(neg_edc_2, edc_inner_r, edc_outer_r, lay_neg_z);
    add_disc_portal(neg_edc_2, edc_inner_r, edc_outer_r, lay_pos_z);

    // Index of first portal
    prev_vol_idx = neg_edc_2.index() - 1;
    next_vol_idx = neg_edc_2.index() + 1;
    layer_pt_index = neg_edc_2.range()[1] - 4;
    surfaces[layer_pt_index].set_edge({beampipe_idx, inv_sf_finder});
    surfaces[++layer_pt_index].set_edge({leaving_world, inv_sf_finder});
    surfaces[++layer_pt_index].set_edge({prev_vol_idx, inv_sf_finder});
    surfaces[++layer_pt_index].set_edge({next_vol_idx, inv_sf_finder});

    //
    // gap layer
    //

    // build gap volume
    gap_neg_z = side * edc_vol_sizes[1].first;
    gap_pos_z = side * edc_vol_sizes[0].second;

    volume_type& neg_edc_gap_2 =
        add_cylinder_volume(edc_inner_r, edc_outer_r, gap_neg_z, gap_pos_z);

    // negative and positive, inner and outer portal surface
    add_cylinder_portal(neg_edc_gap_2, edc_inner_r, gap_neg_z, gap_pos_z);
    add_cylinder_portal(neg_edc_gap_2, edc_outer_r, gap_neg_z, gap_pos_z);
    add_disc_portal(neg_edc_gap_2, edc_inner_r, edc_outer_r, gap_neg_z);
    add_disc_portal(neg_edc_gap_2, edc_inner_r, edc_outer_r, gap_pos_z);

    // Connect portals (first layer, second layer)
    prev_vol_idx = neg_edc_gap_2.index() - 1;
    next_vol_idx = neg_edc_gap_2.index() + 1;
    gap_pt_index = neg_edc_gap_2.range()[1] - 4;
    // cylinder
    surfaces[gap_pt_index].set_edge({beampipe_idx, inv_sf_finder});
    surfaces[++gap_pt_index].set_edge({leaving_world, inv_sf_finder});
    // discs
    surfaces[++gap_pt_index].set_edge({prev_vol_idx, inv_sf_finder});
    surfaces[++gap_pt_index].set_edge({next_vol_idx, inv_sf_finder});

    //
    // third layer
    //
    lay_neg_z = side * edc_vol_sizes[0].second;
    lay_pos_z = side * edc_vol_sizes[0].first;

    volume_type& neg_edc_1 =
        add_cylinder_volume(edc_inner_r, edc_outer_r, lay_neg_z, lay_pos_z);

    // create disc module surfaces
    auto [neg_edc_1_modules, neg_edc_1_trfs, neg_edc_1_masks] = 
                    create_endcap_modules<surface, trapezoid>(
                                                       side * edc_positions[0],
                                                       edc_ring_stagger,
                                                       edc_phi_stagger,
                                                       edc_phi_sub_stagger,
                                                       edc_inner_r, 
                                                       edc_outer_r,
                                                       edc_disc_binning,
                                                       edc_half_y,
                                                       edc_half_x_min_y,
                                                       edc_half_x_max_y,
                                                       edc_tilt);

    // update linking in volumes and surfaces
    update_links(neg_edc_1, neg_edc_1_modules, transforms.size(), trapezoids.size());
    // Append to collections
    surfaces.insert(surfaces.end(), neg_edc_1_modules.begin(), neg_edc_1_modules.end());
    transforms.insert(transforms.end(), neg_edc_1_trfs.begin(), neg_edc_1_trfs.end());
    trapezoids.insert(trapezoids.end(), neg_edc_1_masks.begin(), neg_edc_1_masks.end());

    // negative and positive, inner and outer portal surface
    add_cylinder_portal(neg_edc_1, edc_inner_r, lay_neg_z, lay_pos_z);
    add_cylinder_portal(neg_edc_1, edc_outer_r, lay_neg_z, lay_pos_z);
    add_disc_portal(neg_edc_1, edc_inner_r, edc_outer_r, lay_neg_z);
    add_disc_portal(neg_edc_1, edc_inner_r, edc_outer_r, lay_pos_z);

    // Index of first portal
    prev_vol_idx = neg_edc_1.index() - 1;
    next_vol_idx = neg_edc_1.index() + 1;
    layer_pt_index = neg_edc_1.range()[1] - 4;
    surfaces[layer_pt_index].set_edge({beampipe_idx, inv_sf_finder});
    surfaces[++layer_pt_index].set_edge({leaving_world, inv_sf_finder});
    surfaces[++layer_pt_index].set_edge({prev_vol_idx, inv_sf_finder});
    surfaces[++layer_pt_index].set_edge({next_vol_idx, inv_sf_finder});

    //
    // gap layer
    //

    // build gap volume
    gap_neg_z = side * edc_vol_sizes[0].first;
    gap_pos_z = side * brl_half_z;

    volume_type& neg_edc_gap_1 =
        add_cylinder_volume(edc_inner_r, edc_outer_r, gap_neg_z, gap_pos_z);

    // negative and positive, inner and outer portal surface
    add_cylinder_portal(neg_edc_gap_1, edc_inner_r, gap_neg_z, gap_pos_z);
    add_cylinder_portal(neg_edc_gap_1, edc_outer_r, gap_neg_z, gap_pos_z);
    add_disc_portal(neg_edc_gap_1, edc_inner_r, edc_outer_r, gap_neg_z);
    add_disc_portal(neg_edc_gap_1,  brl_vol_sizes[1].first,  brl_vol_sizes[1].second, gap_pos_z);
    add_disc_portal(neg_edc_gap_1,  brl_vol_sizes[1].second,  brl_vol_sizes[2].first, gap_pos_z);
    add_disc_portal(neg_edc_gap_1,  brl_vol_sizes[2].first,  brl_vol_sizes[2].second, gap_pos_z);
    add_disc_portal(neg_edc_gap_1,  brl_vol_sizes[2].second,  brl_vol_sizes[3].first, gap_pos_z);
    add_disc_portal(neg_edc_gap_1,  brl_vol_sizes[3].first,  brl_vol_sizes[3].second, gap_pos_z);
    add_disc_portal(neg_edc_gap_1,  brl_vol_sizes[3].second,  brl_vol_sizes[4].first, gap_pos_z);
    add_disc_portal(neg_edc_gap_1,  brl_vol_sizes[4].first,  brl_vol_sizes[4].second, gap_pos_z);

    // Connect portals (first layer, second layer)
    prev_vol_idx = neg_edc_gap_1.index() - 1;
    next_vol_idx = neg_edc_gap_1.index() + 1;
    gap_pt_index = neg_edc_gap_1.range()[1] - 10;
    // cylinder
    surfaces[gap_pt_index].set_edge({beampipe_idx, inv_sf_finder});
    surfaces[++gap_pt_index].set_edge({leaving_world, inv_sf_finder});
    // discs
    surfaces[++gap_pt_index].set_edge({prev_vol_idx, inv_sf_finder});
    surfaces[++gap_pt_index].set_edge({next_vol_idx, inv_sf_finder});
    surfaces[++gap_pt_index].set_edge({++next_vol_idx, inv_sf_finder});
    surfaces[++gap_pt_index].set_edge({++next_vol_idx, inv_sf_finder});
    surfaces[++gap_pt_index].set_edge({++next_vol_idx, inv_sf_finder});
    surfaces[++gap_pt_index].set_edge({++next_vol_idx, inv_sf_finder});
    surfaces[++gap_pt_index].set_edge({++next_vol_idx, inv_sf_finder});
    surfaces[++gap_pt_index].set_edge({++next_vol_idx, inv_sf_finder});

    //
    // barrel
    //

    //
    // first layer
    //

    // build the first layer volume
    scalar lay_inner_r = brl_vol_sizes[1].first;
    scalar lay_outer_r = brl_vol_sizes[1].second;

    volume_type& brl_layer_1 =
        add_cylinder_volume(lay_inner_r, lay_outer_r, -brl_half_z, brl_half_z);

    // create module surfaces
    auto [brl_l1_modules, brl_l1_trfs, brl_l1_masks] = create_modules<surface, rectangle>(brl_half_x, brl_half_y, brl_tilt_phi, brl_positions[1], brl_radial_stagger, brl_l_overlap, brl_binning[1]);

    // update linking in volumes and surfaces
    update_links(brl_layer_1, brl_l1_modules, transforms.size(), rectangles.size());
    // Append to collections
    surfaces.insert(surfaces.end(), brl_l1_modules.begin(), brl_l1_modules.end());
    transforms.insert(transforms.end(), brl_l1_trfs.begin(), brl_l1_trfs.end());
    rectangles.insert(rectangles.end(), brl_l1_masks.begin(), brl_l1_masks.end());

    // negative and positive, inner and outer portal surface
    add_cylinder_portal(brl_layer_1, lay_inner_r, -brl_half_z, brl_half_z);
    add_cylinder_portal(brl_layer_1, lay_outer_r, -brl_half_z, brl_half_z);
    add_disc_portal(brl_layer_1, lay_inner_r, lay_outer_r, -brl_half_z);
    add_disc_portal(brl_layer_1, lay_inner_r, lay_outer_r, brl_half_z);

    // Connect portals (beampipe, gap)
    next_vol_idx = brl_layer_1.index() + 1;
    layer_pt_index = brl_layer_1.range()[1] - 4;
    // cylinder
    surfaces[layer_pt_index].set_edge({beampipe_idx, inv_sf_finder});
    surfaces[++layer_pt_index].set_edge({next_vol_idx, inv_sf_finder});
    // discs
    surfaces[++layer_pt_index].set_edge({6, inv_sf_finder});
    surfaces[++layer_pt_index].set_edge({14, inv_sf_finder});

    //
    // gap layer
    //

    // build gap volume (first layer outer r and second layer inner r)
    volume_type& brl_gap_1 = add_cylinder_volume(
         brl_vol_sizes[1].second,  brl_vol_sizes[2].first, -brl_half_z, brl_half_z);

    // negative and positive, inner and outer portal surface
    add_cylinder_portal(brl_gap_1,  brl_vol_sizes[1].second, -brl_half_z, brl_half_z);
    add_cylinder_portal(brl_gap_1, brl_vol_sizes[2].first, -brl_half_z, brl_half_z);
    add_disc_portal(brl_gap_1,  brl_vol_sizes[1].second, brl_vol_sizes[2].first,
                    -brl_half_z);
    add_disc_portal(brl_gap_1,  brl_vol_sizes[1].second, brl_vol_sizes[2].first,
                    brl_half_z);

    // Connect portals (first layer, second layer)
    prev_vol_idx = brl_gap_1.index() - 1;
    next_vol_idx = brl_gap_1.index() + 1;
    gap_pt_index = brl_gap_1.range()[1] - 4;
    // cylinder
    surfaces[gap_pt_index].set_edge({prev_vol_idx, inv_sf_finder});
    surfaces[++gap_pt_index].set_edge({next_vol_idx, inv_sf_finder});
    // discs
    surfaces[++gap_pt_index].set_edge({6, inv_sf_finder});
    surfaces[++gap_pt_index].set_edge({14, inv_sf_finder});

    //
    // second layer
    //

    // build the first layer volume
    lay_inner_r = brl_vol_sizes[2].first;
    lay_outer_r = brl_vol_sizes[2].second;

    volume_type& brl_layer_2 =
        add_cylinder_volume(lay_inner_r, lay_outer_r, -brl_half_z, brl_half_z);

    // create module surfaces
    auto [brl_l2_modules, brl_l2_trfs, brl_l2_masks] = create_modules<surface, rectangle>(brl_half_x, brl_half_y, brl_tilt_phi, brl_positions[2], brl_radial_stagger, brl_l_overlap, brl_binning[2]);

    // update linking in volumes and surfaces
    update_links(brl_layer_2, brl_l2_modules, transforms.size(), rectangles.size());
    // Append to collections
    surfaces.insert(surfaces.end(), brl_l2_modules.begin(), brl_l2_modules.end());
    transforms.insert(transforms.end(), brl_l2_trfs.begin(), brl_l2_trfs.end());
    rectangles.insert(rectangles.end(), brl_l2_masks.begin(), brl_l2_masks.end());

    // negative and positive, inner and outer portal surface
    add_cylinder_portal(brl_layer_2, lay_inner_r, -brl_half_z, brl_half_z);
    add_cylinder_portal(brl_layer_2, lay_outer_r, -brl_half_z, brl_half_z);
    add_disc_portal(brl_layer_2, lay_inner_r, lay_outer_r, -brl_half_z);
    add_disc_portal(brl_layer_2, lay_inner_r, lay_outer_r, brl_half_z);

    // Connect portals (beampipe, gap)
    prev_vol_idx = brl_layer_2.index() - 1;
    next_vol_idx = brl_layer_2.index() + 1;
    layer_pt_index = brl_layer_2.range()[1] - 4;
    // cylinder
    surfaces[layer_pt_index].set_edge({prev_vol_idx, inv_sf_finder});
    surfaces[++layer_pt_index].set_edge({next_vol_idx, inv_sf_finder});
    // discs
    surfaces[++layer_pt_index].set_edge({6, inv_sf_finder});
    surfaces[++layer_pt_index].set_edge({14, inv_sf_finder});

    //
    // gap layer
    //

    // build gap volume (first layer outer r and second layer inner r)
    volume_type& brl_gap_2 = add_cylinder_volume(
         brl_vol_sizes[2].second,  brl_vol_sizes[3].first, -brl_half_z, brl_half_z);

    // negative and positive, inner and outer portal surface
    add_cylinder_portal(brl_gap_2,  brl_vol_sizes[2].second, -brl_half_z, brl_half_z);
    add_cylinder_portal(brl_gap_2, brl_vol_sizes[3].first, -brl_half_z, brl_half_z);
    add_disc_portal(brl_gap_2,  brl_vol_sizes[2].second, brl_vol_sizes[3].first,
                    -brl_half_z);
    add_disc_portal(brl_gap_2,  brl_vol_sizes[2].second, brl_vol_sizes[3].first,
                    brl_half_z);

    // Connect portals (first layer, second layer)
    prev_vol_idx = brl_gap_2.index() - 1;
    next_vol_idx = brl_gap_2.index() + 1;
    gap_pt_index = brl_gap_2.range()[1] - 4;
    // cylinder
    surfaces[gap_pt_index].set_edge({prev_vol_idx, inv_sf_finder});
    surfaces[++gap_pt_index].set_edge({next_vol_idx, inv_sf_finder});
    // discs
    surfaces[++gap_pt_index].set_edge({6, inv_sf_finder});
    surfaces[++gap_pt_index].set_edge({14, inv_sf_finder});

    //
    // third layer
    //

    // build the first layer volume
    lay_inner_r = brl_vol_sizes[3].first;
    lay_outer_r = brl_vol_sizes[3].second;

    volume_type& brl_layer_3 =
        add_cylinder_volume(lay_inner_r, lay_outer_r, -brl_half_z, brl_half_z);

    // create module surfaces
    auto [brl_l3_modules, brl_l3_trfs, brl_l3_masks] = create_modules<surface, rectangle>(brl_half_x, brl_half_y, brl_tilt_phi, brl_positions[3], brl_radial_stagger, brl_l_overlap, brl_binning[3]);

    // update linking in volumes and surfaces
    update_links(brl_layer_3, brl_l3_modules, transforms.size(), rectangles.size());
    // Append to collections
    surfaces.insert(surfaces.end(), brl_l3_modules.begin(), brl_l3_modules.end());
    transforms.insert(transforms.end(), brl_l3_trfs.begin(), brl_l3_trfs.end());
    rectangles.insert(rectangles.end(), brl_l3_masks.begin(), brl_l3_masks.end());

    // negative and positive, inner and outer portal surface
    add_cylinder_portal(brl_layer_3, lay_inner_r, -brl_half_z, brl_half_z);
    add_cylinder_portal(brl_layer_3, lay_outer_r, -brl_half_z, brl_half_z);
    add_disc_portal(brl_layer_3, lay_inner_r, lay_outer_r, -brl_half_z);
    add_disc_portal(brl_layer_3, lay_inner_r, lay_outer_r, brl_half_z);

    // Connect portals (beampipe, gap)
    prev_vol_idx = brl_layer_3.index() - 1;
    next_vol_idx = brl_layer_3.index() + 1;
    layer_pt_index = brl_layer_3.range()[1] - 4;
    // cylinder
    surfaces[layer_pt_index].set_edge({prev_vol_idx, inv_sf_finder});
    surfaces[++layer_pt_index].set_edge({next_vol_idx, inv_sf_finder});
    // discs
    surfaces[++layer_pt_index].set_edge({6, inv_sf_finder});
    surfaces[++layer_pt_index].set_edge({14, inv_sf_finder});

    //
    // gap layer
    //

    // build gap volume (first layer outer r and second layer inner r)
    volume_type& brl_gap_3 = add_cylinder_volume(
         brl_vol_sizes[3].second,  brl_vol_sizes[4].first, -brl_half_z, brl_half_z);

    // negative and positive, inner and outer portal surface
    add_cylinder_portal(brl_gap_3,  brl_vol_sizes[3].second, -brl_half_z, brl_half_z);
    add_cylinder_portal(brl_gap_3, brl_vol_sizes[4].first, -brl_half_z, brl_half_z);
    add_disc_portal(brl_gap_3,  brl_vol_sizes[3].second, brl_vol_sizes[4].first,
                    -brl_half_z);
    add_disc_portal(brl_gap_3,  brl_vol_sizes[3].second, brl_vol_sizes[4].first,
                    brl_half_z);

    // Connect portals (first layer, second layer)
    prev_vol_idx = brl_gap_3.index() - 1;
    next_vol_idx = brl_gap_3.index() + 1;
    gap_pt_index = brl_gap_3.range()[1] - 4;
    // cylinder
    surfaces[gap_pt_index].set_edge({prev_vol_idx, inv_sf_finder});
    surfaces[++gap_pt_index].set_edge({next_vol_idx, inv_sf_finder});
    // discs
    surfaces[++gap_pt_index].set_edge({6, inv_sf_finder});
    surfaces[++gap_pt_index].set_edge({14, inv_sf_finder});

    //
    // fourth layer
    //

    // build the first layer volume
    lay_inner_r = brl_vol_sizes[4].first;
    lay_outer_r = brl_vol_sizes[4].second;

    volume_type& brl_layer_4 =
        add_cylinder_volume(lay_inner_r, lay_outer_r, -brl_half_z, brl_half_z);

    // create module surfaces
    auto [brl_l4_modules, brl_l4_trfs, brl_l4_masks] = create_modules<surface, rectangle>(brl_half_x, brl_half_y, brl_tilt_phi, brl_positions[4], brl_radial_stagger, brl_l_overlap, brl_binning[4]);

    // update linking in volumes and surfaces
    update_links(brl_layer_4, brl_l4_modules, transforms.size(), rectangles.size());
    // Append to collections
    surfaces.insert(surfaces.end(), brl_l4_modules.begin(), brl_l4_modules.end());
    transforms.insert(transforms.end(), brl_l4_trfs.begin(), brl_l4_trfs.end());
    rectangles.insert(rectangles.end(), brl_l4_masks.begin(), brl_l4_masks.end());

    // negative and positive, inner and outer portal surface
    add_cylinder_portal(brl_layer_4, lay_inner_r, -brl_half_z, brl_half_z);
    add_cylinder_portal(brl_layer_4, lay_outer_r, -brl_half_z, brl_half_z);
    add_disc_portal(brl_layer_4, lay_inner_r, lay_outer_r, -brl_half_z);
    add_disc_portal(brl_layer_4, lay_inner_r, lay_outer_r, brl_half_z);

    // Connect portals (beampipe, gap)
    prev_vol_idx = brl_layer_4.index() - 1;
    layer_pt_index = brl_layer_4.range()[1] - 4;
    // cylinder
    surfaces[layer_pt_index].set_edge({prev_vol_idx, inv_sf_finder});
    surfaces[++layer_pt_index].set_edge({leaving_world, inv_sf_finder});
    // discs
    surfaces[++layer_pt_index].set_edge({6, inv_sf_finder});
    surfaces[++layer_pt_index].set_edge({14, inv_sf_finder});

    //
    // positive endcap
    //

    side = 1.;

    //
    // gap layer
    //

    // build gap volume
    gap_neg_z = side * brl_half_z;
    gap_pos_z = side * edc_vol_sizes[0].first;

    volume_type& pos_edc_gap_1 =
        add_cylinder_volume(edc_inner_r, edc_outer_r, gap_neg_z, gap_pos_z);

    // negative and positive, inner and outer portal surface
    add_cylinder_portal(pos_edc_gap_1, edc_inner_r, gap_neg_z, gap_pos_z);
    add_cylinder_portal(pos_edc_gap_1, edc_outer_r, gap_neg_z, gap_pos_z);
    add_disc_portal(pos_edc_gap_1,  brl_vol_sizes[1].first,  brl_vol_sizes[1].second, gap_neg_z);
    add_disc_portal(pos_edc_gap_1,  brl_vol_sizes[1].second,  brl_vol_sizes[2].first, gap_neg_z);
    add_disc_portal(pos_edc_gap_1,  brl_vol_sizes[2].first,  brl_vol_sizes[2].second, gap_neg_z);
    add_disc_portal(pos_edc_gap_1,  brl_vol_sizes[2].second,  brl_vol_sizes[3].first, gap_neg_z);
    add_disc_portal(pos_edc_gap_1,  brl_vol_sizes[3].first,  brl_vol_sizes[3].second, gap_neg_z);
    add_disc_portal(pos_edc_gap_1,  brl_vol_sizes[3].second,  brl_vol_sizes[4].first, gap_neg_z);
    add_disc_portal(pos_edc_gap_1,  brl_vol_sizes[4].first,  brl_vol_sizes[4].second, gap_neg_z);
    add_disc_portal(pos_edc_gap_1, edc_inner_r, edc_outer_r, gap_pos_z);

    // Connect portals (first layer, second layer)
    prev_vol_idx = 7;
    next_vol_idx = pos_edc_gap_1.index() + 1;
    gap_pt_index = pos_edc_gap_1.range()[1] - 10;
    // cylinder
    surfaces[gap_pt_index].set_edge({beampipe_idx, inv_sf_finder});
    surfaces[++gap_pt_index].set_edge({leaving_world, inv_sf_finder});
    // discs
    surfaces[++gap_pt_index].set_edge({prev_vol_idx, inv_sf_finder});
    surfaces[++gap_pt_index].set_edge({++prev_vol_idx, inv_sf_finder});
    surfaces[++gap_pt_index].set_edge({++prev_vol_idx, inv_sf_finder});
    surfaces[++gap_pt_index].set_edge({++prev_vol_idx, inv_sf_finder});
    surfaces[++gap_pt_index].set_edge({++prev_vol_idx, inv_sf_finder});
    surfaces[++gap_pt_index].set_edge({++prev_vol_idx, inv_sf_finder});
    surfaces[++gap_pt_index].set_edge({++prev_vol_idx, inv_sf_finder});
    surfaces[++gap_pt_index].set_edge({next_vol_idx, inv_sf_finder});

    //
    // first layer
    //
    lay_neg_z = side * edc_vol_sizes[0].first;
    lay_pos_z = side * edc_vol_sizes[0].second;

    volume_type& pos_edc_1 =
        add_cylinder_volume(edc_inner_r, edc_outer_r, lay_neg_z, lay_pos_z);

    // create disc module surfaces
    auto [pos_edc_1_modules, pos_edc_1_trfs, pos_edc_1_masks] = 
                    create_endcap_modules<surface, trapezoid>(
                                                       side * edc_positions[0],
                                                       edc_ring_stagger,
                                                       edc_phi_stagger,
                                                       edc_phi_sub_stagger,
                                                       edc_inner_r, 
                                                       edc_outer_r,
                                                       edc_disc_binning,
                                                       edc_half_y,
                                                       edc_half_x_min_y,
                                                       edc_half_x_max_y,
                                                       edc_tilt);

    // update linking in volumes and surfaces
    update_links(pos_edc_1, pos_edc_1_modules, transforms.size(), trapezoids.size());
    // Append to collections
    surfaces.insert(surfaces.end(), pos_edc_1_modules.begin(), pos_edc_1_modules.end());
    transforms.insert(transforms.end(), pos_edc_1_trfs.begin(), pos_edc_1_trfs.end());
    trapezoids.insert(trapezoids.end(), pos_edc_1_masks.begin(), pos_edc_1_masks.end());

    // negative and positive, inner and outer portal surface
    add_cylinder_portal(pos_edc_1, edc_inner_r, lay_neg_z, lay_pos_z);
    add_cylinder_portal(pos_edc_1, edc_outer_r, lay_neg_z, lay_pos_z);
    add_disc_portal(pos_edc_1, edc_inner_r, edc_outer_r, lay_neg_z);
    add_disc_portal(pos_edc_1, edc_inner_r, edc_outer_r, lay_pos_z);

    // Index of first portal
    prev_vol_idx = pos_edc_1.index() - 1;
    next_vol_idx = pos_edc_1.index() + 1;
    layer_pt_index = pos_edc_1.range()[1] - 4;
    surfaces[layer_pt_index].set_edge({beampipe_idx, inv_sf_finder});
    surfaces[++layer_pt_index].set_edge({leaving_world, inv_sf_finder});
    surfaces[++layer_pt_index].set_edge({prev_vol_idx, inv_sf_finder});
    surfaces[++layer_pt_index].set_edge({next_vol_idx, inv_sf_finder});

    //
    // gap layer
    //

    // build gap volume
    gap_neg_z = side * edc_vol_sizes[0].second;
    gap_pos_z = side * edc_vol_sizes[1].first;

    volume_type& pos_edc_gap_2 =
        add_cylinder_volume(edc_inner_r, edc_outer_r, gap_neg_z, gap_pos_z);

    // negative and positive, inner and outer portal surface
    add_cylinder_portal(pos_edc_gap_2, edc_inner_r, gap_neg_z, gap_pos_z);
    add_cylinder_portal(pos_edc_gap_2, edc_outer_r, gap_neg_z, gap_pos_z);
    add_disc_portal(pos_edc_gap_2, edc_inner_r, edc_outer_r, gap_neg_z);
    add_disc_portal(pos_edc_gap_2, edc_inner_r, edc_outer_r, gap_pos_z);

    // Connect portals (first layer, second layer)
    prev_vol_idx = pos_edc_gap_2.index() - 1;
    next_vol_idx = pos_edc_gap_2.index() + 1;
    gap_pt_index = pos_edc_gap_2.range()[1] - 4;
    // cylinder
    surfaces[gap_pt_index].set_edge({beampipe_idx, inv_sf_finder});
    surfaces[++gap_pt_index].set_edge({leaving_world, inv_sf_finder});
    // discs
    surfaces[++gap_pt_index].set_edge({prev_vol_idx, inv_sf_finder});
    surfaces[++gap_pt_index].set_edge({next_vol_idx, inv_sf_finder});

    //
    // second layer
    //
    lay_neg_z = side * edc_vol_sizes[1].first;
    lay_pos_z = side * edc_vol_sizes[1].second;

    volume_type& pos_edc_2 =
        add_cylinder_volume(edc_inner_r, edc_outer_r, lay_neg_z, lay_pos_z);

    // create disc module surfaces
    auto [pos_edc_2_modules, pos_edc_2_trfs, pos_edc_2_masks] = 
                    create_endcap_modules<surface, trapezoid>(
                                                       side * edc_positions[1],
                                                       edc_ring_stagger,
                                                       edc_phi_stagger,
                                                       edc_phi_sub_stagger,
                                                       edc_inner_r, 
                                                       edc_outer_r,
                                                       edc_disc_binning,
                                                       edc_half_y,
                                                       edc_half_x_min_y,
                                                       edc_half_x_max_y,
                                                       edc_tilt);

    // update linking in volumes and surfaces
    update_links(pos_edc_2, pos_edc_2_modules, transforms.size(), trapezoids.size());
    // Append to collections
    surfaces.insert(surfaces.end(), pos_edc_2_modules.begin(), pos_edc_2_modules.end());
    transforms.insert(transforms.end(), pos_edc_2_trfs.begin(), pos_edc_2_trfs.end());
    trapezoids.insert(trapezoids.end(), pos_edc_2_masks.begin(), pos_edc_2_masks.end());

    // negative and positive, inner and outer portal surface
    add_cylinder_portal(pos_edc_2, edc_inner_r, lay_neg_z, lay_pos_z);
    add_cylinder_portal(pos_edc_2, edc_outer_r, lay_neg_z, lay_pos_z);
    add_disc_portal(pos_edc_2, edc_inner_r, edc_outer_r, lay_neg_z);
    add_disc_portal(pos_edc_2, edc_inner_r, edc_outer_r, lay_pos_z);

    // Index of first portal
    prev_vol_idx = pos_edc_2.index() - 1;
    next_vol_idx = pos_edc_2.index() + 1;
    layer_pt_index = pos_edc_2.range()[1] - 4;
    surfaces[layer_pt_index].set_edge({beampipe_idx, inv_sf_finder});
    surfaces[++layer_pt_index].set_edge({leaving_world, inv_sf_finder});
    surfaces[++layer_pt_index].set_edge({prev_vol_idx, inv_sf_finder});
    surfaces[++layer_pt_index].set_edge({next_vol_idx, inv_sf_finder});

    //
    // gap layer
    //

    // build gap volume
    gap_neg_z = side * edc_vol_sizes[1].second;
    gap_pos_z = side * edc_vol_sizes[2].first;

    volume_type& pos_edc_gap_3 =
        add_cylinder_volume(edc_inner_r, edc_outer_r, gap_neg_z, gap_pos_z);

    // negative and positive, inner and outer portal surface
    add_cylinder_portal(pos_edc_gap_3, edc_inner_r,gap_neg_z, gap_pos_z);
    add_cylinder_portal(pos_edc_gap_3, edc_outer_r, gap_neg_z, gap_pos_z);
    add_disc_portal(pos_edc_gap_3, edc_inner_r, edc_outer_r, gap_neg_z);
    add_disc_portal(pos_edc_gap_3, edc_inner_r, edc_outer_r, gap_pos_z);

    // Connect portals (first layer, second layer)
    prev_vol_idx = pos_edc_gap_3.index() - 1;
    next_vol_idx = pos_edc_gap_3.index() + 1;
    gap_pt_index = pos_edc_gap_3.range()[1] - 4;
    // cylinder
    surfaces[gap_pt_index].set_edge({beampipe_idx, inv_sf_finder});
    surfaces[++gap_pt_index].set_edge({leaving_world, inv_sf_finder});
    // discs
    surfaces[++gap_pt_index].set_edge({prev_vol_idx, inv_sf_finder});
    surfaces[++gap_pt_index].set_edge({next_vol_idx, inv_sf_finder});

    //
    // third layer
    //

    // build the outermost volume
    lay_neg_z = side * edc_vol_sizes[2].first;
    lay_pos_z = side * edc_vol_sizes[2].second;

    volume_type& pos_edc_3 =
        add_cylinder_volume(edc_inner_r, edc_outer_r, lay_neg_z, lay_pos_z);

    // create disc module surfaces
    auto [pos_edc_3_modules, pos_edc_3_trfs, pos_edc_3_masks] = 
                    create_endcap_modules<surface, trapezoid>(
                                                       side * edc_positions[2],
                                                       edc_ring_stagger,
                                                       edc_phi_stagger,
                                                       edc_phi_sub_stagger,
                                                       edc_inner_r, 
                                                       edc_outer_r,
                                                       edc_disc_binning,
                                                       edc_half_y,
                                                       edc_half_x_min_y,
                                                       edc_half_x_max_y,
                                                       edc_tilt);

    // update linking in volumes and surfaces
    update_links(pos_edc_3, pos_edc_3_modules, transforms.size(), trapezoids.size());
    // Append to collections
    surfaces.insert(surfaces.end(), pos_edc_3_modules.begin(), pos_edc_3_modules.end());
    transforms.insert(transforms.end(), pos_edc_3_trfs.begin(), pos_edc_3_trfs.end());
    trapezoids.insert(trapezoids.end(), pos_edc_3_masks.begin(), pos_edc_3_masks.end());

    // negative and positive, inner and outer portal surface
    add_cylinder_portal(pos_edc_3, edc_inner_r, lay_neg_z, lay_pos_z);
    add_cylinder_portal(pos_edc_3, edc_outer_r, lay_neg_z, lay_pos_z);
    add_disc_portal(pos_edc_3, edc_inner_r, edc_outer_r, lay_neg_z);
    add_disc_portal(pos_edc_3, edc_inner_r, edc_outer_r, lay_pos_z);

    // Index of first portal
    prev_vol_idx = pos_edc_3.index() - 1;
    layer_pt_index = pos_edc_3.range()[1] - 4;
    surfaces[layer_pt_index].set_edge({beampipe_idx, inv_sf_finder});
    surfaces[++layer_pt_index].set_edge({leaving_world, inv_sf_finder});
    surfaces[++layer_pt_index].set_edge({prev_vol_idx, inv_sf_finder});
    surfaces[++layer_pt_index].set_edge({leaving_world, inv_sf_finder});

    // Return all geometry containers
    using geometry_t = toy_geometry<volume_type, surface>;

    geometry_t geo(resource);
    geo.add_volumes(std::move(volumes));
    geo.add_objects(std::move(surfaces));

    // First, put data into the detector interface
    /*mask_store<dtuple, dvector, discs, cylinders, rectangles> masks;
    // populate mask store
    masks.add_masks(discs);
    masks.add_masks(cylinders);
    masks.add_masks(rectangles);*/
    auto masks = std::make_tuple<rectangle_container, trapezoid_container, annulus_container, cylinder_container, disc_container>(
        std::move(rectangles), std::move(trapezoids), std::move(annuli),
        std::move(cylinders), std::move(discs));

    auto d = toy_detector<decltype(geo), decltype(masks),
                          decltype(transforms)::value_type, toy_grid>(
        std::move(geo), std::move(masks), std::move(transforms));

    return std::move(d);
}

}  // namespace detray
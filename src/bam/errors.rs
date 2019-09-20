use snafu::Snafu;

pub type Result<T, E = Error> = std::result::Result<T, E>;

#[derive(Snafu, Debug)]
#[snafu(visibility = "pub")]
pub enum Error {
    // migrate all BAM related errors here, see bcf::errors for a blueprint.
}
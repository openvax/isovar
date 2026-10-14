# Complete report ledgers as release assets (#501)

The completed forward summary exceeds GitHub's 100 MiB Git-file limit. Its
witnesses also appear in the complete event ledger. Preserve that ledger's
schema and every value; distribute the oversized file as a checksummed,
versioned release asset. Keep the original full report in its verified portable
archive and prepare a separate repository copy with updated Markdown links.

Validate every source report file against its existing checksum manifest before
preparation. Copy the ledger with bounded buffers and verify its complete bytes.
Publish a download inventory recording its logical filename, URL, size and SHA256,
plus the source report manifest checksum. Generate a new checksum manifest for
the repository copy, including the externally hosted ledger's logical filename.
Never overwrite different assets or publish an incomplete output directory.

Test exact asset bytes, source immutability, Markdown links, local/external
checksum resolution, changed-source rejection, existing-output rejection and
the release-asset size limit. After the standard clean-master PyPI deployment,
upload the prepared asset to the matching GitHub release and verify the downloaded
bytes before claiming publication complete. The running recovery validators and
queued archive/report helpers remain frozen; this distribution step follows them.

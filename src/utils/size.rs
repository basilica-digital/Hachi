//! Human-readable byte-size parsing/formatting for CLI arguments.
//!
//! Accepts inputs like `"1024"`, `"1024B"`, `"10KB"`, `"1MiB"`, `"1gb"`,
//! `"1Tb"` (case-insensitive; an optional `B`/`iB`/`ib` suffix and
//! intervening whitespace are permitted). All unit multipliers are
//! powers of 1024 — `1 KB == 1 KiB == 2^10`, etc. — because the
//! surrounding cryptographic scheme only cares about powers of two.

/// Parse a human-readable size string into a byte count.
pub fn parse_size(s: &str) -> Result<u64, String> {
    let s = s.trim();
    if s.is_empty() {
        return Err("empty size string".into());
    }

    // Split into leading numeric part and trailing unit part.
    let unit_start = s
        .char_indices()
        .find(|(_, c)| c.is_ascii_alphabetic())
        .map(|(i, _)| i)
        .unwrap_or(s.len());
    let (num_str, unit_str) = s.split_at(unit_start);

    let num: u64 = num_str
        .trim()
        .parse()
        .map_err(|e| format!("invalid number '{}': {e}", num_str.trim()))?;

    let mult: u64 = match unit_str.trim().to_ascii_lowercase().as_str() {
        "" | "b" => 1,
        "k" | "kb" | "kib" => 1 << 10,
        "m" | "mb" | "mib" => 1 << 20,
        "g" | "gb" | "gib" => 1 << 30,
        "t" | "tb" | "tib" => 1 << 40,
        other => {
            return Err(format!(
                "unknown size unit '{other}' (use B/KB/MB/GB/TB, case-insensitive)"
            ))
        }
    };

    num.checked_mul(mult)
        .ok_or_else(|| format!("size overflow: {num} * {mult}"))
}

/// Format a byte count using the largest binary unit that divides it evenly.
pub fn format_size(bytes: u64) -> String {
    const KIB: u64 = 1 << 10;
    const MIB: u64 = 1 << 20;
    const GIB: u64 = 1 << 30;
    const TIB: u64 = 1 << 40;

    if bytes >= TIB && bytes % TIB == 0 {
        format!("{} TiB", bytes / TIB)
    } else if bytes >= GIB && bytes % GIB == 0 {
        format!("{} GiB", bytes / GIB)
    } else if bytes >= MIB && bytes % MIB == 0 {
        format!("{} MiB", bytes / MIB)
    } else if bytes >= KIB && bytes % KIB == 0 {
        format!("{} KiB", bytes / KIB)
    } else {
        format!("{bytes} B")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn plain_bytes() {
        assert_eq!(parse_size("1024").unwrap(), 1024);
        assert_eq!(parse_size("1024B").unwrap(), 1024);
        assert_eq!(parse_size("1024b").unwrap(), 1024);
        assert_eq!(parse_size("0").unwrap(), 0);
    }

    #[test]
    fn kilobytes() {
        assert_eq!(parse_size("10KB").unwrap(), 10 << 10);
        assert_eq!(parse_size("10kb").unwrap(), 10 << 10);
        assert_eq!(parse_size("10KiB").unwrap(), 10 << 10);
        assert_eq!(parse_size("10 KiB").unwrap(), 10 << 10);
        assert_eq!(parse_size("1k").unwrap(), 1 << 10);
    }

    #[test]
    fn megabytes() {
        assert_eq!(parse_size("1MB").unwrap(), 1 << 20);
        assert_eq!(parse_size("1mb").unwrap(), 1 << 20);
        assert_eq!(parse_size("1MiB").unwrap(), 1 << 20);
    }

    #[test]
    fn gigabytes() {
        assert_eq!(parse_size("1GB").unwrap(), 1 << 30);
        assert_eq!(parse_size("10Gb").unwrap(), 10u64 << 30);
        assert_eq!(parse_size("4GiB").unwrap(), 4u64 << 30);
    }

    #[test]
    fn terabytes() {
        assert_eq!(parse_size("1TB").unwrap(), 1u64 << 40);
        assert_eq!(parse_size("1Tb").unwrap(), 1u64 << 40);
        assert_eq!(parse_size("1TiB").unwrap(), 1u64 << 40);
    }

    #[test]
    fn rejects_bad_input() {
        assert!(parse_size("").is_err());
        assert!(parse_size("   ").is_err());
        assert!(parse_size("10ZB").is_err());
        assert!(parse_size("1.5GB").is_err());
        assert!(parse_size("-10KB").is_err());
        assert!(parse_size("abc").is_err());
    }

    #[test]
    fn round_trip() {
        for s in ["1024", "10 KiB", "4GiB", "1TB"] {
            let bytes = parse_size(s).unwrap();
            let formatted = format_size(bytes);
            assert_eq!(parse_size(&formatted).unwrap(), bytes);
        }
    }
}

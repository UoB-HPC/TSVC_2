#!/usr/bin/env bash
set -euo pipefail

apt_options=(-o Acquire::Retries=5 -o Acquire::http::Timeout=30 -o Acquire::https::Timeout=30)
curl_options=(--fail --silent --show-error --location --retry 5 --retry-all-errors --retry-delay 30 --connect-timeout 30)

apt-get "${apt_options[@]}" update
apt-get "${apt_options[@]}" install -y --no-install-recommends make curl ca-certificates gpg g++

if [[ -n ${CI_APT_SOURCE:-} ]]; then
  printf 'Downloading repository key: %s\n' "$CI_APT_KEY"
  curl "${curl_options[@]}" "$CI_APT_KEY" |
    gpg --dearmor --yes -o /usr/share/keyrings/tsvc-toolchain.gpg
  printf 'Adding package repository: %s\n' "$CI_APT_SOURCE"
  printf 'deb [signed-by=/usr/share/keyrings/tsvc-toolchain.gpg] %s\n' \
    "$CI_APT_SOURCE" > /etc/apt/sources.list.d/tsvc-toolchain.list
  repo_host=${CI_APT_SOURCE#*://}
  repo_host=${repo_host%%/*}
  printf 'Package: *\nPin: origin %s\nPin-Priority: 600\n' "$repo_host" \
    > /etc/apt/preferences.d/tsvc-toolchain
  apt-get "${apt_options[@]}" update
fi

read -r -a packages <<< "$CI_PACKAGES"
if ((${#packages[@]})); then
  apt-get "${apt_options[@]}" install -y --no-install-recommends "${packages[@]}"
fi

read -r -a cc_command <<< "$CC"
read -r -a cxx_command <<< "$CXX"
"${cc_command[0]}" --version
"${cxx_command[0]}" --version
make --version

#!/bin/sh

set -e

[ -z "${GH_TOKEN}" ] && exit 0
[ "${TRAVIS_BRANCH}" != "master" ] && exit 0

git clone -b gh-pages https://${GH_TOKEN}@github.com/${TRAVIS_REPO_SLUG}.git book-output
cd book-output
cp -r ../_book/* ./
git add --all *
git commit -m"Update the book"
git push origin gh-pages

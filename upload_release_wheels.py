import requests
import zipfile
import io
import os
import argparse

min_py3 = 7
max_py3 = 10


def download_file(url, folder_name):
    r = requests.get(url)
    z = zipfile.ZipFile(io.BytesIO(r.content))
    z.extractall(folder_name)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--upload", required=True, type=int, help="0: only download wheels"
                                                              "1: download and upload wheels on Pypi")
    args = vars(ap.parse_args())
    upload = args["upload"]

    wheels_folder = "./tmp_folder"
    for py3_ver in range(min_py3, max_py3 + 1):
        print("Downloading Windows wheels for python 3." + str(py3_ver))
        wheel_url = "https://ci.appveyor.com/api/projects/stfbnc/fathon/artifacts/wheels.zip?" \
                    "job=Environment:%20CIBW_BUILD=cp3" + str(py3_ver) + "-win_amd64"
        download_file(wheel_url, wheels_folder)

    print("Downloading latest artifacts from github")
    os.system("gh run download --dir " + wheels_folder)

    for f in os.listdir(wheels_folder):
        if os.path.isdir(os.path.join(wheels_folder, f)):
            curr_dir = os.path.join(wheels_folder, f)
            os.system("mv " + curr_dir + "/* " + wheels_folder)
            os.system("rm -rf " + curr_dir)

    if upload == 1:
        print("Uploading on Pypi")
        os.system("twine upload --verbose --skip-existing " + wheels_folder + "/*")
        os.system("rm " + wheels_folder + "/*")

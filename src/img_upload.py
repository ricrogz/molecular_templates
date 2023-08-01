import json
import os
from base64 import b64encode

from requests_oauthlib import OAuth2Session

ALBUM_HASH = 'Eo6pEWB'
IMGUR_TOKEN_REFRESH_URL = 'https://api.imgur.com/oauth2/token'
IMGUR_UPLOAD_URL = "https://api.imgur.com/3/upload.json"
UPLOAD_RETRIES = 5

IMGUR_CLIENT_ID = os.environ.get('IMGUR_CLIENT_ID', None)
IMGUR_CLIENT_SECRET = os.environ.get('IMGUR_CLIENT_SECRET', None)
IMGUR_TOKEN = json.loads(os.environ.get('IMGUR_TOKEN', ''))


def export_imgur_token(token):
    IMGUR_TOKEN = token
    if gh_output := os.environ.get('GITHUB_OUTPUT', ''):
        with open(gh_output, 'a') as f:
            f.write(f'imgur_token={json.dumps(IMGUR_TOKEN)}')


def upload_img(fpath, title):
    session = OAuth2Session(client_id=IMGUR_CLIENT_ID,
                            token=IMGUR_TOKEN,
                            auto_refresh_url=IMGUR_TOKEN_REFRESH_URL,
                            auto_refresh_kwargs={
                                'client_id': IMGUR_CLIENT_ID,
                                'client_secret': IMGUR_CLIENT_SECRET
                            },
                            token_updater=print)

    img = b64encode(open(fpath, 'rb').read())
    fname = os.path.basename(fpath)
    data = {'image': img, 'type': 'base64', 'name': fname, 'title': title}
    if ALBUM_HASH:
        data['album'] = ALBUM_HASH

    for _ in range(UPLOAD_RETRIES):
        request = session.post(IMGUR_UPLOAD_URL, data=data)
        if request.ok:
            break
    if not request.ok:
        raise RuntimeError(request.text)

    return request.json()['data']['link']
